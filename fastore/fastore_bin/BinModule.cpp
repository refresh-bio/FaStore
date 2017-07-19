/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "Globals.h"

#include <vector>
#include <array>
#include <iostream>

#include "BinModule.h"
#include "FastqStream.h"
#include "FastqParser.h"
#include "FastqPacker.h"
#include "FastqCategorizer.h"
#include "BinFile.h"
#include "BinOperator.h"
#include "Exception.h"
#include "Thread.h"


void BinModuleSE::Fastq2Bin(const std::vector<std::string> &inFastqFiles_, const std::string &outBinFile_,
							const BinModuleConfig& config_, uint32 threadNum_,
							bool compressedInput_, bool verboseMode_)
{
	// TODO: try/catch to free resources
	//
	IFastqStreamReaderSE* fastqFile = NULL;
	if (compressedInput_)
		fastqFile = new MultiFastqFileReaderGzSE(inFastqFiles_);
	else
		fastqFile = new MultiFastqFileReaderSE(inFastqFiles_);


	BinFileWriter binFile;
	binFile.StartCompress(outBinFile_, config_);

	const uint32 minimizersCount = config_.minimizer.TotalMinimizersCount();
	if (threadNum_ > 1)
	{
		typedef TFastqChunkReader<IFastqStreamReaderSE, FastqChunkCollectionSE> FastqChunkReader;
		typedef typename FastqChunkReader::FastqChunkPool FastqChunkPool;
		typedef typename FastqChunkReader::FastqChunkQueue FastqChunkQueue;

		FastqChunkPool* fastqPool = NULL;
		FastqChunkQueue* fastqQueue = NULL;
		BinaryPartsPool* binPool = NULL;
		BinaryPartsQueue* binQueue = NULL;

		FastqChunkReader* fastqReader = NULL;
		BinChunkWriter* binWriter = NULL;

		const uint32 partNum = threadNum_ + (threadNum_ >> 2);
		fastqPool = new FastqChunkPool(partNum, config_.fastqBlockSize);
		fastqQueue = new FastqChunkQueue(partNum, 1);

		binPool = new BinaryPartsPool(partNum, minimizersCount);
		binQueue = new BinaryPartsQueue(partNum, threadNum_);

		fastqReader = new FastqChunkReader(fastqFile, fastqQueue, fastqPool);
		binWriter = new BinChunkWriter(&binFile, binQueue, binPool);

		// launch stuff
		//
		mt::thread readerThread(mt::ref(*fastqReader));

		std::vector<IOperator*> operators;
		operators.resize(threadNum_);

		std::vector<mt::thread> opThreadGroup;

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			operators[i] = new BinEncoderSE(config_, fastqQueue, fastqPool, binQueue, binPool);
			opThreadGroup.push_back(mt::thread(mt::ref(*operators[i])));
		}

		(*binWriter)();

		readerThread.join();

		for (mt::thread& t : opThreadGroup)
		{
			t.join();
		}

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			delete operators[i];
		}

		TFREE(binWriter);
		TFREE(fastqReader);

		TFREE(binQueue);
		TFREE(binPool);
		TFREE(fastqQueue);
		TFREE(fastqPool);
	}
	else
	{
		FastqRecordsParserSE parser(config_.archiveType.readsHaveHeaders);
		FastqCategorizerSE categorizer(config_.minimizer, config_.minFilter, config_.catParams);
		FastqRecordsPackerSE packer(config_);

		FastqChunkCollectionSE fastqChunk(config_.fastqBlockSize);
		std::vector<FastqRecord> reads;
		reads.resize(1 << 10);

		std::map<uint32, FastqRecordsPtrBin> dnaBins;
		BinaryBinBlock binBins;

#if (DEV_DEBUG_MODE)
		std::vector<FastqRecord> outReads;
		FastqChunkCollectionSE outChunk;
#endif

		FastqRawBlockStats parseStats;
		while (fastqFile->ReadNextChunk(fastqChunk))
		{
			parseStats.Clear();
			parser.ParseFrom(fastqChunk, reads, parseStats, config_.headParams.preserveComments);
			ASSERT(reads.size() > 0);

			categorizer.Categorize(reads, dnaBins);

			binBins.Clear();
			packer.PackToBins(dnaBins, binBins);

#if (DEV_DEBUG_MODE && 0)
			packer.UnpackFromBin(binBins, outReads, outChunk);
			ASSERT(outReads.size() == reads.size());

			// we need to sort here
			FastqComparator comparator;
			std::sort(reads.begin(), reads.end(), comparator);
			std::sort(outReads.begin(), outReads.end(), comparator);

			std::array<char, FastqRecord::MaxSeqLen> rcBuffer;
			FastqRecord rcRec;
			rcRec.seq = rcBuffer.data();

			for (uint64 i = 0; i < reads.size(); ++i)
			{
				const FastqRecord &r1 = reads[i];
				FastqRecord r2 = outReads[i];

				if (r1.IsReadReverse())
				{
					// WARN: use CopyTo/CopyFrom
					r2.ComputeRC(rcRec);
					std::copy(rcRec.seq, rcRec.seq + rcRec.seqLen, r2.seq);
				}

				ASSERT(std::equal(r1.seq, r1.seq + r1.seqLen, r2.seq));
			}
#endif

			// update stats
			//
			binBins.stats.Update(parseStats);
			binFile.WriteNextBlock(&binBins);
		}
	}

	binFile.FinishCompress();

	if (verboseMode_)
	{
#if 0
		IBinFile::BinFileFooter bf = binFile.GetFileFooter();

		std::cout << "#Signatures count: " << bf.binOffsets.size()
				  << " / " <<  bf.params.minimizer.TotalMinimizersCount() << " + 1 " << std::endl;

		std::cout << "#Records distribution in bins by signature:\n";
		std::cout << "signature_id\tsignature\trecords_count\td_meta_min\td_meta_max\td_dna_min\td_dna_max\n";
		for (std::map<uint32, IBinFile::BinFileFooter::BinInfo>::const_iterator iB = bf.binOffsets.begin();
			 iB != bf.binOffsets.end();
			 iB++)
		{
			char sig[64];
			sig[bf.params.minimizer.signatureLen] = 0;
			bf.params.minimizer.GenerateMinimizer(iB->first, sig);

			std::cout << iB->first << '\t'
					  << sig << '\t'
					  << iB->second.totalRecordsCount << '\t';
			/*
			if (iB->second.blocksMetaData.size() > 1)
			{
				std::cout << iB->second.minMetaDeltaValue << '\t'
						  << iB->second.maxMetaDeltaValue << '\t'
						  << iB->second.minDnaDeltaValue << '\t'
						  << iB->second.maxDnaDeltaValue << '\n';
			}
			else
			{
				std::cout << "*\t*\t*\t*\n";
			}
			*/
		}

#endif

		// TODO: here we can print some extra quality stats verbose info
		//
		//

		const auto& header = binFile.GetFileHeader();
		const auto& footer = binFile.GetFileFooter();

		uint64 totalRawDna = 0, totalRawId = 0;
		uint64 totalPackedDna = 0, totalPackedQua = 0, totalPackedId = 0;

		for (const auto& off : footer.binOffsets)
		{
			totalRawDna += off.second.totalRawDnaSize;
			totalRawId += off.second.totalRawHeadSize;
			totalPackedDna += off.second.totalDnaSize;
			totalPackedQua += off.second.totalQuaSize;
			totalPackedId += off.second.totalHeadSize;
		}
		std::cout << "DNA: " << totalRawDna << " --> " << totalPackedDna << std::endl;
		std::cout << "QUA: " << totalRawDna << " --> " << totalPackedQua << std::endl;
		std::cout << "ID: " << totalRawId << " --> " << totalPackedId << std::endl;

		std::cout << "Records count: " << header.recordsCount << std::endl;
		std::cout << "File footer size: " << header.footerSize << std::endl;

	}

	delete fastqFile;
}


void BinModuleSE::Bin2Dna(const std::string &inBinFile_, const std::string &outFile_)
{
	// TODO: try/catch to free resources
	//
	ASSERT(outFile_.size() != 0);

	BinFileReader binFile;

	BinModuleConfig config;
	binFile.StartDecompress(inBinFile_, config);

	FastqFileWriterSE dnaFile(outFile_);
	FastqChunkCollectionSE fastqChunk(config.fastqBlockSize >> 1);			// WARNING! --- here can be a BUG
	FastqChunkCollectionSE outChunk(config.fastqBlockSize >> 1);			// WARNING! --- here can be a BUG
	FastqRecordsPackerSE packer(config);
	FastqRecordsParserSE parser(config.archiveType.readsHaveHeaders);

	BinaryBinBlock binBin;
	std::vector<FastqRecord> reads;

	while (binFile.ReadNextBlock(&binBin))
	{
		packer.UnpackFromBin(binBin, reads, *fastqChunk.chunks[0], false);
		parser.ParseTo(reads, outChunk);

		dnaFile.WriteNextChunk(outChunk);
	}

	dnaFile.Close();
	binFile.FinishDecompress();
}

void BinModulePE::Fastq2Bin(const std::vector<std::string>& inFastqFiles_1_,
							const std::vector<std::string>& inFastqFiles_2_,
							const std::string & outBinFile_, const BinModuleConfig& config_,
							uint32 threadNum_, bool compressedInput_, bool verboseMode_)
{

	// TODO: try/catch to free resources
	//
	ASSERT(!inFastqFiles_1_.empty());
	ASSERT(inFastqFiles_1_.size() == inFastqFiles_2_.size());

	IFastqStreamReaderPE* fastqFile = NULL;
	if (compressedInput_)
		fastqFile = new MultiFastqFileReaderGzPE(inFastqFiles_1_, inFastqFiles_2_);
	else
		fastqFile = new MultiFastqFileReaderPE(inFastqFiles_1_, inFastqFiles_2_);

	BinFileWriter binFile;
	binFile.StartCompress(outBinFile_, config_);

	const uint32 minimizersCount = config_.minimizer.TotalMinimizersCount();

	if (threadNum_ > 1)
	{
		typedef TFastqChunkReader<IFastqStreamReaderPE, FastqChunkCollectionPE> FastqChunkReader;
		typedef typename FastqChunkReader::FastqChunkPool FastqChunkPool;
		typedef typename FastqChunkReader::FastqChunkQueue FastqChunkQueue;

		FastqChunkPool* fastqPool = NULL;
		FastqChunkQueue* fastqQueue = NULL;
		BinaryPartsPool* binPool = NULL;
		BinaryPartsQueue* binQueue = NULL;

		FastqChunkReader* fastqReader = NULL;
		BinChunkWriter* binWriter = NULL;

		const uint32 partNum = threadNum_ * 2;
		fastqPool = new FastqChunkPool(partNum, config_.fastqBlockSize);
		fastqQueue = new FastqChunkQueue(partNum, 1);

		binPool = new BinaryPartsPool(partNum, minimizersCount);
		binQueue = new BinaryPartsQueue(partNum, threadNum_);

		fastqReader = new FastqChunkReader(fastqFile, fastqQueue, fastqPool);
		binWriter = new BinChunkWriter(&binFile, binQueue, binPool);

		// launch stuff
		//
		mt::thread readerThread(mt::ref(*fastqReader));

		std::vector<IOperator*> operators;
		operators.resize(threadNum_);

		std::vector<mt::thread> opThreadGroup;

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			operators[i] = new BinEncoderPE(config_, fastqQueue, fastqPool, binQueue, binPool);
			opThreadGroup.push_back(mt::thread(mt::ref(*operators[i])));
		}

		(*binWriter)();

		readerThread.join();

		for (mt::thread& t : opThreadGroup)
		{
			t.join();
		}

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			delete operators[i];
		}

		TFREE(binWriter);
		TFREE(fastqReader);

		TFREE(binQueue);
		TFREE(binPool);
		TFREE(fastqQueue);
		TFREE(fastqPool);
	}
	else
	{

		FastqRecordsParserPE parser(config_.archiveType.readsHaveHeaders, config_.headParams.preserveComments);
		FastqCategorizerPE categorizer(config_.minimizer, config_.minFilter, config_.catParams);
		FastqRecordsPackerPE packer(config_);

		FastqChunkCollectionPE inputChunk(config_.fastqBlockSize);
		FastqChunk fastqChunk(config_.fastqBlockSize);
		std::vector<FastqRecord> records;
		records.resize(1 << 10);

		std::map<uint32, FastqRecordsPtrBin> dnaBins;
		BinaryBinBlock binBins;

		FastqRawBlockStats stats;
		while (fastqFile->ReadNextChunk(inputChunk))		// it just extracts RAW FASTQ file chunks
		{
			stats.Clear();
			parser.ParseFrom(inputChunk, records, stats, config_.headParams.preserveComments);

			dnaBins.clear();
			categorizer.Categorize(records, dnaBins);

			binBins.Clear();
			packer.PackToBins(dnaBins, binBins);

			// update stats
			//
			binBins.stats.Update(stats);
			binFile.WriteNextBlock(&binBins);
		}
	}

	binFile.FinishCompress();


	if (verboseMode_)
	{
		const auto& header = binFile.GetFileHeader();
		const auto& footer = binFile.GetFileFooter();

		uint64 totalRawDna = 0, totalRawId = 0;
		uint64 totalPackedDna = 0, totalPackedQua = 0, totalPackedId = 0;

		for (const auto& off : footer.binOffsets)
		{
			totalRawDna += off.second.totalRawDnaSize;
			totalRawId += off.second.totalRawHeadSize;
			totalPackedDna += off.second.totalDnaSize;
			totalPackedQua += off.second.totalQuaSize;
			totalPackedId += off.second.totalHeadSize;
		}
		std::cout << "DNA: " << totalRawDna << " --> " << totalPackedDna << std::endl;
		std::cout << "QUA: " << totalRawDna << " --> " << totalPackedQua << std::endl;
		std::cout << "ID: " << totalRawId << " --> " << totalPackedId << std::endl;

		std::cout << "Records count: " << header.recordsCount << std::endl;
		std::cout << "File footer size: " << header.footerSize << std::endl;
	}

	/*
	if (verboseMode_)
	{
		std::vector<uint64> recordCounts;
		binFile.GetBinStats(recordCounts);

		std::cout << "Signatures count: " << recordCounts.size() << std::endl;
		std::cout << "Records distribution in bins by signature:\n";
		for (uint32 i = 0; i < recordCounts.size(); ++i)
		{
			if (recordCounts[i] > 0)
				std::cout << i << " : " << recordCounts[i] << '\n';
		}
		std::cout << std::endl;
	}
	*/

	delete fastqFile;
}



void BinModulePE::Bin2Dna(const std::string & inBinFile_,
						  const std::string & outFile1_,
						  const std::string & outFile2_)
{
	// TODO: try/catch to free resources
	//
	BinModuleConfig config;
	BinFileReader binFile;
	binFile.StartDecompress(inBinFile_, config);

	const auto& headInfo = binFile.GetFileFooter().headData;

	FastqFileWriterPE fqFile(outFile1_, outFile2_);
	FastqChunk outChunk(config.fastqBlockSize >> 1);
	FastqChunkCollectionPE fastqChunk(config.fastqBlockSize >> 1);

	FastqRecordsPackerPE packer(config);
	FastqRecordsParserPE parser(config.archiveType.readsHaveHeaders, headInfo.pairedEndFieldIdx);

	BinaryBinBlock binBin;
	std::vector<FastqRecord> reads;

	while (binFile.ReadNextBlock(&binBin))
	{
		packer.UnpackFromBin(binBin, reads, outChunk, false);
		parser.ParseTo(reads, fastqChunk);

		fqFile.WriteNextChunk(fastqChunk);
	}

	fqFile.Close();
	binFile.FinishDecompress();
}

