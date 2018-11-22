/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "Globals.h"

#include <vector>
#include <map>

#include "BinOperator.h"
#include "FastqStream.h"
#include "BinFile.h"
#include "FastqParser.h"
#include "FastqCategorizer.h"
#include "FastqPacker.h"

#include <iostream>


void BinChunkWriter::Run()
{
	int64 partId = 0;
	BinaryBinBlock* part = NULL;

	uint64 partsProcessed = 0;
	while (partsQueue->Pop(partId, part))
	{
		// here PartId is not important
		partsStream->WriteNextBlock(part);

		// release the memory
		//
		part->Reset();
		partsPool->Release(part);
		part = NULL;

		if (verboseMode)
		{
			partsProcessed++;

			std::cerr << '\r' << "Parts processed: " << partsProcessed;

			if (totalPartsCount > 0)
				std::cerr << " (" << partsProcessed * 100 / totalPartsCount << "%)";
			std::cerr << std::flush;
		}
	}
}


struct BinBuffer
{
	static const uint32 MinRecordsToStore = 64;

	std::vector<FastqRecord> records;
	DataChunk buffer;

	BinBuffer(uint32 maxRecords_, uint64 bufferSize_)
		:	buffer(bufferSize_)
	{
		records.reserve(maxRecords_);
	}
};


void BinEncoderSE::Run()
{
	std::map<uint32, BinBuffer*> binBuffers;

	FastqRecordsPackerSE packer(binConfig);
	FastqCategorizerSE categorizer(binConfig.minimizer, binConfig.minFilter, binConfig.catParams);

	int64 partId = 0;
	FastqChunkCollectionSE* fqPart = NULL;
	BinaryBinBlock* binPart = NULL;

	std::map<uint32, FastqRecordsPtrBin> dnaBins;
	std::vector<FastqRecord> reads;
	reads.resize(1 << 10);

	FastqRecordsParserSE parser(binConfig.archiveType.readsHaveHeaders);

	FastqRawBlockStats stats;
	while (fqPartsQueue->Pop(partId, fqPart))			// different types
	{
		// TIP: when processing small files, stats need to be cleared at the end of each bin processing, 
		// as the stats will be lost if the bin will be empty after post-processing
		//stats.Clear();
		parser.ParseFrom(*fqPart, reads, stats, binConfig.headParams.preserveComments);				// different types
		ASSERT(!reads.empty());

		categorizer.Categorize(reads, dnaBins);


		// TODO: move this to some common code-base
		//


		// check the bins
		//
		std::vector<uint32> binsToClear;
		for (auto i = dnaBins.begin(); i != dnaBins.end(); )
		{
			if (i->second.records.size() == 0)
			{
				i = dnaBins.erase(i);
				continue;
			}

			if (i->second.records.size() < BinBuffer::MinRecordsToStore)
			{
				const uint32 sig = i->first;

				auto& records = i->second.records;

				// try to merge if
				//
				if (binBuffers.count(sig) != 0 && binBuffers.at(sig)->records.size() + records.size() >= BinBuffer::MinRecordsToStore)
				{
					auto& bufRecords = binBuffers.at(sig)->records;

					for (auto& rec : bufRecords)
						records.push_back((FastqRecord*)&rec);

					// erase the buffer bin
					//
					binsToClear.push_back(sig);


					// continue the iteration
					//
					i++;
				}
				else
				{
					const FastqRecord& r0 = *records.back();	// back to have the max id number
					const bool usesQuality = r0.qua != NULL;
					const bool usesHeaders = r0.head != NULL;
					const uint64 approxReadSize = (uint64)r0.seqLen * (1 + (uint64)usesQuality) + (uint64)usesHeaders*r0.headLen*1.2;
					const uint64 approxChunkSize = records.size() * approxReadSize;

					if (binBuffers.count(sig) == 0)
						binBuffers[sig] = new BinBuffer(BinBuffer::MinRecordsToStore, MAX(BinBuffer::MinRecordsToStore * approxReadSize, approxChunkSize));

					DataChunk& buffer = binBuffers.at(sig)->buffer;
					auto& bufRecords = binBuffers.at(sig)->records;

					ASSERT(bufRecords.size() + records.size() < BinBuffer::MinRecordsToStore);

					// copy data
					//
					for (const auto& rec : records)
					{
						const uint64 readSize = (uint64)rec->seqLen * (1 + (uint64)usesQuality) + (uint64)rec->headLen;
						ASSERT(buffer.size + readSize < buffer.data.Size());

						char* bufferPtr = (char*)buffer.data.Pointer() + buffer.size;

						bufRecords.push_back(FastqRecord(*rec));
						FastqRecord& bufRec = bufRecords.back();

						bufRec.seq = bufferPtr;
						if (usesQuality)
						{
							bufRec.qua = bufferPtr + bufRec.seqLen;
						}

						if (usesHeaders)
						{
							bufRec.head = bufferPtr + (bufRec.seqLen * (1 + (uint32)usesQuality));
						}

						bufRec.CopyFrom(*rec, true);

						buffer.size += readSize;
					}


					// erase the bin
					//
					i = dnaBins.erase(i);
				}

				// WARN: remember to kill read count checkup at the end
			}
			else
			{
				i++;
			}
		}


		// do not process empty bins (they were not empty, but now are after post processing)
		if (dnaBins.empty())
			continue;

		partId++;
		binPartsPool->Acquire(binPart);

		packer.PackToBins(dnaBins, *binPart);

		// set stats of the processed part
		// WARN: the stats are 'approximated' as we're also performing bins filtering
		//
		// TODO: swap()
		binPart->stats.Clear();
		binPart->stats.Update(stats);

		fqPartsPool->Release(fqPart);					// this one uses different type

		ASSERT(binPart->descriptors.size() > 0);
		binPartsQueue->Push(partId, binPart);


		// TODO: move to some common code-base
		//


		// cleanup the buffers
		//
		if (binsToClear.size() > 0)
		{
			for (uint32 sig : binsToClear)
			{
				BinBuffer* buf = binBuffers.at(sig);

#if EXTRA_MEM_OPT
				delete buf;
				binBuffers.erase(sig);
#else
				buf->records.clear();
				buf->buffer.size = 0;
#endif
			}
		}


		// TODO: reuse this one
		//
		dnaBins.clear();

		stats.Clear();
	}

	// TODO: move to some common code-base
	//

	// store the remaining reads in the buffer
	//
	// WARN: do not update stats here, as we have already them calulcated while parsing
	//
	{
		dnaBins.clear();
		const uint32 nBinId = binConfig.minimizer.TotalMinimizersCount();

		FastqRecordBuffer rcRec;

		for (auto iBin = binBuffers.begin(); iBin != binBuffers.end(); ++iBin)
		{
			if (iBin->second->records.size() == 0)
				continue;

			// filter the small bins
			//
			if (iBin->second->records.size() < binConfig.catParams.minBlockBinSize)
			{
				auto& dnaBin = dnaBins[nBinId];
				auto& records = dnaBin.records;

				for (FastqRecord& rec : iBin->second->records)
				{
					rec.minimPos = 0;
					if (rec.IsReadReverse())
					{
						rec.ComputeRC(rcRec);			// TODO: optimize join
						rec.CopyFrom(rcRec);
						rec.SetReadReverse(false);
					}
					records.push_back(&rec);
				}

				// here we collect stats
				//
				dnaBin.stats.minSeqLen = dnaBin.stats.maxSeqLen = records.front()->seqLen;
				ASSERT(dnaBin.stats.minSeqLen > 0);
			}
			else
			{
#if DEV_DEBUG_MODE
				char minString[64];
				binConfig.minimizer.GenerateMinimizer(iBin->first, minString);
#endif

				auto& dnaBin = dnaBins[iBin->first];
				auto& records = dnaBin.records;

				for (auto& rec : iBin->second->records)
				{
#if DEV_DEBUG_MODE
					const char* minSeq = std::search(rec.seq, rec.seq + rec.seqLen,
													   minString, minString + binConfig.minimizer.signatureLen);

					ASSERT(minSeq != rec.seq + rec.seqLen || iBin->first == nBinId);
					ASSERT(minSeq == rec.seq + rec.minimPos || iBin->first == nBinId);
#endif
					records.push_back(&rec);
				}

				dnaBin.stats.minSeqLen = dnaBin.stats.maxSeqLen = records.front()->seqLen;
				ASSERT(dnaBin.stats.minSeqLen > 0);
			}
		}
	}


	if (!dnaBins.empty())
	{
		partId = 0;
		binPartsPool->Acquire(binPart);

		packer.PackToBins(dnaBins, *binPart);

		binPart->stats.Update(stats);

		ASSERT(binPart->descriptors.size() > 0);
		binPartsQueue->Push(partId, binPart);
	}
	binPartsQueue->SetCompleted();


	// cleanup the buffer
	//
	for (auto iBin : binBuffers)
	{
		delete iBin.second;
	}
}


void BinEncoderPE::Run()
{
	std::map<uint32, BinBuffer*> binBuffers;

	// those ones need to be templatized
	//
	FastqRecordsPackerPE packer(binConfig);
	FastqCategorizerPE categorizer(binConfig.minimizer, binConfig.minFilter, binConfig.catParams);
	FastqChunkCollectionPE* fqPart = NULL;
	FastqRecordsParserPE parser(binConfig.archiveType.readsHaveHeaders);
	//
	// //

	int64 partId = 0;
	BinaryBinBlock* binPart = NULL;

	std::map<uint32, FastqRecordsPtrBin> dnaBins;
	std::vector<FastqRecord> reads;
	reads.resize(1 << 10);

	FastqRawBlockStats stats;
	while (fqPartsQueue->Pop(partId, fqPart))
	{
		// TODO: templatize + add parser proxy to use one code base
		//

		stats.Clear();
		parser.ParseFrom(*fqPart, reads, stats, binConfig.headParams.preserveComments);
		ASSERT(!reads.empty());

		categorizer.Categorize(reads, dnaBins);

		// check the bins
		//
		std::vector<uint32> binsToClear;
		for (auto i = dnaBins.begin(); i != dnaBins.end(); )
		{
			if (i->second.records.size() == 0)
			{
				i = dnaBins.erase(i);
				continue;
			}

			if (i->second.records.size() < BinBuffer::MinRecordsToStore)
			{
				const uint32 sig = i->first;

				auto& records = i->second.records;

				// try to merge if
				//
				if (binBuffers.count(sig) != 0 && binBuffers.at(sig)->records.size() + records.size() >= BinBuffer::MinRecordsToStore)
				{
					auto& bufRecords = binBuffers.at(sig)->records;

					for (auto& rec : bufRecords)
						records.push_back((FastqRecord*)&rec);

					// erase the buffer bin
					//
					binsToClear.push_back(sig);


					// continue the iteration
					//
					i++;
				}
				else
				{
					const FastqRecord& r0 = *records.back();	// back to have the max id number
					const bool usesQuality = r0.qua != NULL;
					const bool usesHeaders = r0.head != NULL;
					const uint64 approxRecordSize = (uint64)(r0.seqLen + r0.auxLen) * (1 + (uint64)usesQuality) + (uint64)usesHeaders * (uint64)r0.headLen * 1.2;
					const uint64 approxChunkSize = records.size() * approxRecordSize;


					if (binBuffers.count(sig) == 0)
						binBuffers[sig] = new BinBuffer(BinBuffer::MinRecordsToStore, BinBuffer::MinRecordsToStore * approxRecordSize);

					DataChunk& buffer = binBuffers.at(sig)->buffer;
					auto& bufRecords = binBuffers.at(sig)->records;

					ASSERT(bufRecords.size() + records.size() < BinBuffer::MinRecordsToStore);
					ASSERT(buffer.size + approxChunkSize < BinBuffer::MinRecordsToStore * approxRecordSize);


					// copy data
					//
					for (const auto& rec : records)
					{
						const uint64 readSize = (uint64)(rec->seqLen + rec->auxLen) * (1 + (uint64)usesQuality) + (uint64)rec->headLen;
						char* bufferPtr = (char*)buffer.data.Pointer() + buffer.size;

						bufRecords.push_back(FastqRecord(*rec));
						FastqRecord& bufRec = bufRecords.back();

						bufRec.seq = bufferPtr;
						if (usesQuality)
						{
							bufRec.qua = bufferPtr + bufRec.seqLen + bufRec.auxLen;
						}
						if (usesHeaders)
						{
							bufRec.head = bufferPtr + (bufRec.seqLen + bufRec.auxLen) * (1 + (uint32)usesQuality);
						}

						bufRec.CopyFrom(*rec, true);

						buffer.size += readSize;
					}


					// erase the bin
					//
					i = dnaBins.erase(i);
				}

				// WARN: remember to kill read count checkup at the end
			}
			else
			{
				i++;
			}
		}

		partId++;
		binPartsPool->Acquire(binPart);
		packer.PackToBins(dnaBins, *binPart);

		fqPartsPool->Release(fqPart);

		// update stats
		//
		binPart->stats.Clear();
		binPart->stats.Update(stats);

		binPartsQueue->Push(partId, binPart);


		// cleanup the buffers
		//
		if (binsToClear.size() > 0)
		{
			for (uint32 sig : binsToClear)
			{
				BinBuffer* buf = binBuffers.at(sig);
#if EXTRA_MEM_OPT
				delete buf;
				binBuffers.erase(sig);
#else
				buf->records.clear();
				buf->buffer.size = 0;
#endif
			}
		}


		// clear the bins -- TODO: reuse
		//
		dnaBins.clear();
	}


	// store the remaining reads in the buffer
	//
	// WARN: do not update the stats as they were already computed
	// while parsing
	//
	{
		dnaBins.clear();

		FastqRecordBuffer rcRec;
		for (auto iBin = binBuffers.begin(); iBin != binBuffers.end(); ++iBin)
		{
			if (iBin->second->records.size() == 0)
				continue;

			// filter the small bins
			//
			if (iBin->second->records.size() < binConfig.catParams.minBlockBinSize)
			{
				auto& dnaBin = dnaBins[binConfig.minimizer.SignatureN()];
				auto& records = dnaBin.records;

				for (FastqRecord& rec : iBin->second->records)
				{
					rec.minimPos = 0;

					if (rec.IsReadReverse())
					{
						rec.ComputeRC(rcRec);
						rec.CopyFrom(rcRec);					// TODO: optmimize join
						rec.SetReadReverse(false);
					}

					if (rec.IsPairSwapped())
						rec.SwapReads();

					records.push_back(&rec);
				}

				// here we collect stats
				//
				dnaBin.stats.minSeqLen = dnaBin.stats.maxSeqLen = records.front()->seqLen;
				ASSERT(dnaBin.stats.minSeqLen > 0);
			}
			else
			{
#if DEV_DEBUG_MODE
				std::array<char, MAX_SIGNATURE_LEN> minString;
				binConfig.minimizer.GenerateMinimizer(iBin->first, minString.data());
#endif
				auto& dnaBin = dnaBins[iBin->first];
				auto& records = dnaBin.records;

				for (auto& rec : iBin->second->records)
				{
#if DEV_DEBUG_MODE
					const char* minSeq = std::search(rec.seq,
													   rec.seq + rec.seqLen,
													   minString.data(),
													   minString.data() + binConfig.minimizer.signatureLen);
					ASSERT(minSeq != rec.seq + rec.seqLen || iBin->first == binConfig.minimizer.SignatureN());
					ASSERT(minSeq == rec.seq + rec.minimPos || iBin->first == binConfig.minimizer.SignatureN());
#endif

					records.push_back(&rec);
				}

				dnaBin.stats.minSeqLen = dnaBin.stats.maxSeqLen = records.front()->seqLen;
				ASSERT(dnaBin.stats.minSeqLen > 0);
			}
		}
	}


	if (!dnaBins.empty())
	{
		partId = 0;

		binPartsPool->Acquire(binPart);
		packer.PackToBins(dnaBins, *binPart);

		binPartsQueue->Push(partId, binPart);
	}

	binPartsQueue->SetCompleted();


	// cleanup the buffer
	//
	for (auto iBin : binBuffers)
	{
		delete iBin.second;
	}
}
