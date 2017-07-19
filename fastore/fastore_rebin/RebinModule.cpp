/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "../fastore_bin/Globals.h"

#include <vector>
#include <iostream>

#include "RebinModule.h"
#include "RebinOperator.h"
#include "DnaRebalancer.h"
#include "NodesPacker.h"

#include "../fastore_bin/FastqCategorizer.h"

#include "../fastore_bin/BinFile.h"
#include "../fastore_bin/BinOperator.h"
#include "../fastore_bin/Exception.h"
#include "../fastore_bin/Thread.h"

#include "../fastore_pack/BinFileExtractor.h"



#include <stdio.h>

void RebinModule::Bin2Bin(const std::string &inBinFile_,
							const std::string &outBinFile_,
							const BinBalanceParameters& params_,
							uint32 threadsNum_,
							bool verboseMode_)
{
	BinModuleConfig conf;
	BinFileExtractor* extractor = new BinFileExtractor(params_.minBinSizeToExtract);		// here uses default minimum bin size

	extractor->StartDecompress(inBinFile_, conf);
	const bool pairedEnd = conf.archiveType.readType == ArchiveType::READ_PE;
	conf.binningLevel = int_log(params_.signatureParity, 2);

	BinFileWriter* writer = new BinFileWriter();
	writer->StartCompress(outBinFile_, conf);

	const QualityCompressionData& quaCompData = extractor->GetFileFooter().quaData;
	const auto& headData = extractor->GetFileFooter().headData;


	// count bin frequencies
	//
	BinBalanceParameters params(params_);
	const auto descriptors = extractor->GetBlockDescriptors(true);

	const uint32 nBinId = conf.minimizer.TotalMinimizersCount();
	const uint32 totalBinsCount = descriptors.size();

	if (params.minBinSizeToCategorize > 0)
	{
		params.validBinSignatures.resize(conf.minimizer.TotalMinimizersCount(), false);

		for (auto iDesc : descriptors)
		{
			if (iDesc.second->totalRecordsCount >= params.minBinSizeToCategorize)
				params.validBinSignatures[iDesc.first] = true;
		}
	}
	else
	{
		params.validBinSignatures.resize(conf.minimizer.TotalMinimizersCount(), true);
	}

	// TODO: optimization, here we should also applpy the valid signatures mask
	//
	// .. .. .. ..

	std::unique_ptr<IFastqNodesPacker> packer(!pairedEnd
		? (IFastqNodesPacker*)(new FastqNodesPackerSE(conf))
		: (IFastqNodesPacker*)(new FastqNodesPackerPE(conf)));

	std::unique_ptr<FastqRecordsPackerSE> rawPacker(!pairedEnd
		? new FastqRecordsPackerSE(conf)
		: new FastqRecordsPackerPE(conf));

	DnaRebalancer rebalancer(conf.minimizer, params, pairedEnd);

	BinaryBinBlock binBin;
	RebinWorkBuffer binBuffer;

	FastqRecordBinStats stats;

	// the small bins and N just keep intact
	//
#if DEV_DEBUG_MODE
	if (verboseMode_)
		std::cerr << "Processing small bins" << std::endl;
#endif

	FastqRecordBuffer rcRec;

	while (extractor->ExtractNextSmallBin(binBin))
	{
		if (binBin.auxDescriptors.size() > 1)
		{
			// TODO: optimize merger: do not unpack, just merge raw data or
			// change header of the block to being a multi-part bin
			binBuffer.Reset();
			packer->UnpackFromBin(binBin, binBuffer.reads,
								 *binBuffer.rebinCtx.graph,
								 stats, binBuffer.dataBuffer, false);

			// reverse reads
			//
			for (FastqRecord& rec : binBuffer.reads)
			{
				if (rec.IsReadReverse())
				{
					rec.ComputeRC(rcRec);
					rec.CopyFrom(rcRec);
					rec.SetReadReverse(false);
				}

				if (rec.IsPairSwapped())
					rec.SwapReads();

				rec.minimPos = 0;
			}

			rawPacker->PackToBin(binBuffer.reads, binBin, nBinId);
		}

		writer->WriteNextBlock(&binBin);
	}


	// TODO: also try to re-categorize Ns
	//
#if DEV_DEBUG_MODE
	if (verboseMode_)
		std::cerr << "Processing N bin" << std::endl;
#endif

	if (extractor->ExtractNBin(binBin))
	{
		if (binBin.auxDescriptors.size() > 1)
		{
			// TODO: optimize merger: do not unpack, just merge raw data or
			// change header of the block to being a multi-part bin
			binBuffer.Reset();
			packer->UnpackFromBin(binBin, binBuffer.reads,
								 *binBuffer.rebinCtx.graph,
								 stats, binBuffer.dataBuffer, false);

			// reverse reads
			//
			for (FastqRecord& rec : binBuffer.reads)
			{
				if (rec.IsReadReverse())
				{
					rec.ComputeRC(rcRec);
					rec.CopyFrom(rcRec);
					rec.SetReadReverse(false);
				}

				if (rec.IsPairSwapped())
					rec.SwapReads();

				rec.minimPos = 0;
			}

			rawPacker->PackToBin(binBuffer.reads, binBin, nBinId);
		}

		writer->WriteNextBlock(&binBin);
	}

#if DEV_DEBUG_MODE
	if (verboseMode_)
		std::cerr << "Processing standard bins" << std::endl;
#endif

	if (threadsNum_ > 1)
	{
		// TODO: UPDATE ME
		//
		const uint32 partNum = threadsNum_ + (threadsNum_ >> 2);
		const uint64 inBufferSize = 1 << 8;
		const uint64 outBufferSize = 1 << 8;

		MinimizerPartsPool* inPool = new MinimizerPartsPool(partNum, inBufferSize);
		MinimizerPartsQueue* inQueue = new MinimizerPartsQueue(partNum, 1);

		BinaryPartsPool* outPool = new BinaryPartsPool(partNum, outBufferSize);
		BinaryPartsQueue* outQueue = new BinaryPartsQueue(partNum, threadsNum_);

		BinPartsExtractor* inReader = new BinPartsExtractor(extractor, inQueue, inPool);
		BinChunkWriter* outWriter = new BinChunkWriter(writer,
													   outQueue,
													   outPool,
													   verboseMode_,
													   totalBinsCount);


		// launch threads
		//
		mt::thread readerThread(mt::ref(*inReader));

		std::vector<IOperator*> operators;

#ifdef USE_BOOST_THREAD
		boost::thread_group opThreadGroup;

		for (uint32 i = 0; i < threadsNum_; ++i)
		{
			IOperator* op = new BinBalancer(conf.minimizer, params, conf.quaParams,
										   inQueue, inPool,
										   outQueue, outPool, pairedEnd);
			operators.push_back(op);
			opThreadGroup.create_thread(mt::ref(*operators[i]));
		}

		(*outWriter)();

		readerThread.join();
		opThreadGroup.join_all();


#else
		std::vector<mt::thread> opThreadGroup;

		for (uint32 i = 0; i < threadsNum_; ++i)
		{
			IOperator* op = new BinBalancer(conf, params, inQueue, inPool, outQueue, outPool);
			operators.push_back(op);
			opThreadGroup.push_back(mt::thread(mt::ref(*op)));
		}

		(*outWriter)();

		readerThread.join();

		for (mt::thread& t : opThreadGroup)
			t.join();


#endif

		for (IOperator* op : operators)
			delete op;

		TFREE(outWriter);
		TFREE(inReader);

		TFREE(outQueue);
		TFREE(outPool);
		TFREE(inQueue);
		TFREE(inPool);
	}
	else
	{
		uint32 processedBins = 0;
		while (extractor->ExtractNextStdBin(binBin))
		{
			ASSERT(binBin.metaSize > 0);
			ASSERT(binBin.signature != 0);

			const uint32 signatureId = binBin.signature;

#if (DEV_DEBUG_MODE)
			if (verboseMode_)
				std::cout << signatureId << " " << std::flush;
#endif

			if (BinBalanceParameters::IsSignatureValid(signatureId, params_.signatureParity))
			{
				binBuffer.Reset();
				packer->UnpackFromBin(binBin, binBuffer.reads,
									 *binBuffer.rebinCtx.graph, stats,
									 binBuffer.dataBuffer, false);
				ASSERT(binBin.rawDnaSize > 0);

				// here we can already release binbin
				//
				const uint64 inRawReadsCount = binBuffer.reads.size();

				rebalancer.Rebalance(binBuffer.rebinCtx, binBuffer.nodesMap, signatureId);

				packer->PackToBins(binBuffer.nodesMap, binBin);

				uint64 outRawReadsCount = 0;
				for (const auto& desc : binBin.descriptors)
					outRawReadsCount += desc.second.recordsCount;
				ASSERT(outRawReadsCount == inRawReadsCount);

#if (DEV_DEBUG_MODE)
			if (verboseMode_)
				std::cout << "+\n" << std::flush;
#endif
			}
			else
			{
				// TODO: optimize merger: do not unpack, to save memory: just merge raw data or
				// change header of the block to being a multi-part bin
				if (binBin.auxDescriptors.size() > 1)
				{
					binBuffer.Reset();
					packer->UnpackFromBin(binBin, binBuffer.reads,
										 *binBuffer.rebinCtx.graph, stats,
										 binBuffer.dataBuffer, false);
					const uint64 inRawReadsCount = binBuffer.reads.size();

					packer->PackToBin(*binBuffer.rebinCtx.graph, binBin, signatureId);

					uint64 outRawReadsCount = 0;
					for (const auto& desc : binBin.auxDescriptors)
						outRawReadsCount += desc.recordsCount;
					ASSERT(outRawReadsCount == inRawReadsCount);
				}

#if (DEV_DEBUG_MODE)
				if (verboseMode_)
					std::cout << "-\n" << std::flush;
#endif
			}

			writer->WriteNextBlock(&binBin);

			processedBins++;

			if (verboseMode_)
			{
				std::cerr << '\r' << signatureId << " : " << processedBins * 100 / totalBinsCount << "%" << std::flush;
			}
		}
	}

	extractor->FinishDecompress();

	// store the initial FASTQ stats (discarding the calculated while processing)
	//
	writer->SetQualityCompressionData(quaCompData);			// WARN: be careful here about potential memory leak / corruption
	writer->SetHeaderCompressionData(headData);
	writer->FinishCompress();

	if (verboseMode_)
	{
		const IBinFile::BinFileFooter& outBf = writer->GetFileFooter();
		const IBinFile::BinFileFooter& inBf = extractor->GetFileFooter();

		std::cout << "#Signatures count: " << outBf.binOffsets.size()
				  << " ( " << inBf.binOffsets.size() << " ) "
				  << " / " <<  outBf.params.minimizer.TotalMinimizersCount() << " + 1 " << std::endl;

		std::cout << "#Records distribution in bins by signature:\n";
		std::cout << "signature\trecords_count\trecords_diff\td_meta_min\td_meta_max\td_dna_min\td_dna_max\n";
		for (const auto& iB : inBf.binOffsets)
		{
			std::array<char, FastqRecord::MaxSeqLen> sig;
			sig[outBf.params.minimizer.signatureLen] = 0;
			outBf.params.minimizer.GenerateMinimizer(iB.first, sig.data());

			std::cout << iB.first << '\t'
					  << sig.data() << '\t';

			if (outBf.binOffsets.count(iB.first) == 0)
			{
				// the bin has been eliminated
				std::cout << "\t*\t*\t*\t*\t*\t*\n";
				continue;
			}

			const IBinFile::BinFileFooter::BinInfo& oB = outBf.binOffsets.at(iB.first);

			int64 rdiff = oB.totalRecordsCount - iB.second.totalRecordsCount;

			std::cout << oB.totalRecordsCount << '\t'
					  << rdiff << '\t';
			/*
			if (oB.blocksMetaData.size() > 1)
			{
				std::cout << oB.minMetaDeltaValue << '\t'
						  << oB.maxMetaDeltaValue << '\t'
						  << oB.minDnaDeltaValue << '\t'
						  << oB.maxDnaDeltaValue << '\n';
			}
			else
			{
				std::cout << "*\t*\t*\t*\n";
			}
			*/
		}
		std::cout << std::endl;
	}

	if (verboseMode_)
	{
		std::cerr << "\n" << std::flush;
	}

	delete writer;

	// HACK!!!
	QvzCodebook* q = (QvzCodebook *)(&extractor->GetFileFooter().quaData.codebook);
	q->qlist = NULL;

	delete extractor;
}


void RebinModule::Bin2Dna(const std::string &inBinFile_,
						  const std::vector<std::string> &outFiles_)
{
	BinFileReader binFile;
	BinModuleConfig config;
	binFile.StartDecompress(inBinFile_, config);

	FastqChunk inChunk(config.fastqBlockSize >> 1);			// WARNING! --- here can be a BUG

	std::unique_ptr<IFastqStreamWriterSE> writer;
	std::unique_ptr<IFastqChunkCollection> outChunk;
	std::unique_ptr<IFastqNodesPacker> packer;
	std::unique_ptr<IRecordsParser> parser;

	if (config.archiveType.readType == ArchiveType::READ_SE)
	{
		writer.reset(new FastqFileWriterSE(outFiles_[0]));
		outChunk.reset(new FastqChunkCollectionSE(config.fastqBlockSize >> 1));
		packer.reset(new FastqNodesPackerSE(config));
		parser.reset(new FastqRecordsParserSE(config.archiveType.readsHaveHeaders));
	}
	else
	{
		ASSERT(outFiles_.size() == 2);
		writer.reset(new FastqFileWriterPE(outFiles_[0], outFiles_[1]));
		outChunk.reset(new FastqChunkCollectionPE(config.fastqBlockSize >> 1));
		packer.reset(new FastqNodesPackerPE(config));
		parser.reset(new FastqRecordsParserPE(config.archiveType.readsHaveHeaders,
											  binFile.GetFileFooter().headData.pairedEndFieldIdx));
	}



	BinaryBinBlock binBin;
	std::vector<FastqRecord> reads;
	GraphEncodingContext graph;
	FastqRecordBinStats stats;

	//uint64 totalBinsCount = binFile.GetFileFooter().binOffsets.size();
	//uint64 binNum = 0;
	while (binFile.ReadNextBlock(&binBin))
	{
		// clear the buffers
		//
		reads.clear();
		graph.Clear();
		inChunk.Reset();

#if EXTRA_MEM_OPT
		reads.shrink_to_fit();
#endif

		// unpack the reads and parse to output chunk
		//
		packer->UnpackFromBin(binBin, reads, graph, stats, inChunk, false);
		parser->ParseTo(reads, *outChunk);

		// write the reads
		//
		writer->WriteNextChunk(*outChunk);

		//if (verboseMode_)
		//std::cerr << '\r' << binBin.signature << " : " << ++binNum * 100 / totalBinsCount << "%" << std::flush;
	}

	writer->Close();
	binFile.FinishDecompress();
}

