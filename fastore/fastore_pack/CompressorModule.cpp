/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "../core/Globals.h"

#include <iostream>
#include <algorithm>

#include "CompressorModule.h"
#include "ArchiveFile.h"
#include "FastqCompressor.h"
#include "CompressorOperator.h"
#include "Params.h"

#include "../core/FastqPacker.h"
#include "../core/FastqParser.h"
#include "../core/Thread.h"
#include "../core/NodesPacker.h"
#include "../core/FastqStream.h"
#include "../core/FastqParser.h"
#include "../qvz/QVZ.h"
#include "../fastore_rebin/BinFileExtractor.h"
#include "../fastore_rebin/Params.h"


void CompressorModuleSE::Bin2Dnarch(const std::string &inBinFile_, const std::string &outArchiveFile_,
								const CompressorParams& compParams_, const CompressorAuxParams& auxParams_,
								uint32 threadsNum_, bool verboseMode_)
{
	BinModuleConfig binConf;
	BinFileExtractor* extractor = new BinFileExtractor(compParams_.extractor.minBinSize);

	extractor->StartDecompress(inBinFile_, binConf);
	ASSERT(binConf.archiveType.readType != ArchiveType::READ_PE);

	ArchiveFileWriter::ArchiveConfig archConf;
	archConf.archType = binConf.archiveType;
	archConf.minParams = binConf.minimizer;
	archConf.quaParams = binConf.quaParams;

	ArchiveFileWriter* dnarch = new ArchiveFileWriter();
	dnarch->StartCompress(outArchiveFile_, archConf);

	const uint32 totalBinsCount = extractor->GetBlockDescriptors(true).size();
	uint64 nRecordsCount = 0;


	// local compression params used for other sub-modules
	//
	CompressorParams params(compParams_);
	params.archType = archConf.archType;
	params.minimizer = binConf.minimizer;
	params.quality = binConf.quaParams;

	// here we can already calculate global QVZ codebooks
	//
	const QualityCompressionData& globalQuaData = extractor->GetFileFooter().quaData;		// be careful with shared pointers
	// ...

	const auto& headData = extractor->GetFileFooter().headData;


#if (DEV_DEBUG_MODE)
	if (verboseMode_) std::cout << "Processing small bins...\n" << std::flush;
#endif

	CompressedFastqBlockStats stats;
	{	
		FastqCompressor compressor(params, globalQuaData, headData, auxParams_);
		FastqNodesPackerSE packer(binConf);

		// TODO: integrate 3 of those
		std::vector<FastqRecord> reads;
		PackContext mainPackCtx;
		FastqWorkBuffersSE workBuffers;

		CompressedFastqBlock compBin;
		BinaryBinBlock binBin;

		uint64 totalDnaBufferSize = 0;
		uint64 totalHeadBufferSize = 0;
		auto descriptors = extractor->GetBlockDescriptors(false);
		for (const auto& desc : descriptors)
		{
			totalDnaBufferSize += desc.second->totalRawDnaSize;
			totalHeadBufferSize += desc.second->totalRawHeadSize;
			nRecordsCount += desc.second->totalRecordsCount;
		}

		auto nBinDesc = extractor->GetNBlockDescriptor();
		if (nBinDesc.second != NULL)
		{
			totalDnaBufferSize += nBinDesc.second->totalRawDnaSize;
			totalHeadBufferSize += nBinDesc.second->totalRawHeadSize;
			nRecordsCount += nBinDesc.second->totalRecordsCount;
		}

		const uint64 bufferPreallocSize = totalDnaBufferSize * 2
				+ (uint64)binConf.archiveType.readsHaveHeaders * totalHeadBufferSize;
		if (workBuffers.fastqBuffer.data.Size() < bufferPreallocSize)		// support for quality
			workBuffers.fastqBuffer.data.Extend(bufferPreallocSize);

		// WARN: we need to reserve the memory for reads as
		// later we'll be using pointers to records while depacking nodes
		if (reads.size() < nRecordsCount)
			reads.reserve(nRecordsCount);


		// extract and unpack small bins
		//
		//IFastqChunkCollection tmpChunk;

		while (extractor->ExtractNextSmallBin(binBin))
		{
			ASSERT(binBin.metaSize != 0);

			packer.UnpackFromBin(binBin, reads, *mainPackCtx.graph, mainPackCtx.stats, workBuffers.fastqBuffer, true);
		}

		if (extractor->ExtractNBin(binBin))
		{
			packer.UnpackFromBin(binBin, reads, *mainPackCtx.graph, mainPackCtx.stats, workBuffers.fastqBuffer, true);
		}


		// un-reverse-compliment records
		//
		if (reads.size() > 0)
		{
			FastqRecordBuffer rcRec;

			for (FastqRecord& rec : reads)
			{
				if (rec.IsReadReverse())
				{
					rec.ComputeRC(rcRec);
					rec.CopyFrom(rcRec);
					rec.SetReadReverse(false);
				}
				rec.minimPos = 0;
			}


			// check
			//for (auto & rec : reads)
			//{
			//	ASSERT(rec.seq[0] == 'A' || rec.seq[0] == 'C' || rec.seq[0] == 'G' || rec.seq[0] == 'T' || rec.seq[0] == 'N');
			//}

			// compress all bins together
			//
			const uint32 nSignature = params.minimizer.SignatureN();
			compressor.Compress(reads, mainPackCtx, nSignature,
								   totalDnaBufferSize, workBuffers.fastqWorkBin,
								   compBin);

			dnarch->WriteNextBin(compBin.dataBuffer, nSignature);

			stats = compBin.stats;
		}
	}

	if (threadsNum_ > 1)
	{
		const uint32 partNum = threadsNum_ + (threadsNum_ >> 2);
		const uint64 dnaBufferSize = 1 << 8;
		const uint64 outBufferSize = 1 << 8;


		MinimizerPartsPool* inPool = new MinimizerPartsPool(partNum, dnaBufferSize);
		MinimizerPartsQueue* inQueue = new MinimizerPartsQueue(partNum, 1);

		CompressedFastqBlockPool* outPool = new CompressedFastqBlockPool(partNum, outBufferSize);
		CompressedFastqBlockQueue* outQueue = new CompressedFastqBlockQueue(partNum, threadsNum_);

		BinPartsExtractor* inReader = new BinPartsExtractor(extractor, inQueue, inPool);

		ArchivePartsWriter* outWriter = new ArchivePartsWriter(dnarch,
															 outQueue,
															 outPool,
															 verboseMode_,
															 totalBinsCount);

		// update stats from preprocessing
		//
		outWriter->GetStats().Update(stats);


		// TODO : reuse the first part of the buffered data from pool
		//


		// launch stuff
		//
		mt::thread readerThread(mt::ref(*inReader));

		std::vector<IOperator*> operators;
		operators.resize(threadsNum_);

		std::vector<mt::thread> opThreadGroup;

		for (IOperator*& op : operators)
		{
			op = new BinPartsCompressor(params, auxParams_, binConf, globalQuaData, headData,
										inQueue, inPool, outQueue, outPool);
			opThreadGroup.push_back(mt::thread(mt::ref(*op)));
		}

		(*outWriter)();

		readerThread.join();

		for (mt::thread& t : opThreadGroup)
		{
			t.join();
		}

		stats = outWriter->GetStats();

		for (uint32 i = 0; i < threadsNum_; ++i)
		{
			delete operators[i];
		}

		TFREE(outWriter);
		TFREE(inReader);

		TFREE(outQueue);
		TFREE(outPool);
		TFREE(inQueue);
		TFREE(inPool);
	}
	else
	{
		FastqCompressor compressor(params, globalQuaData, headData, auxParams_);
		FastqNodesPackerDynSE packer(binConf);

		CompressedFastqBlock compBin;
		FastqWorkBuffersSE workBuffers;
		std::vector<FastqRecord> reads;
		PackContext mainPackCtx;

		IFastqChunkCollection tmpChunk;

		BinaryBinBlock binBin;

#if(DEV_DEBUG_MODE)
		FastqCompressor decompressor(params, globalQuaData);
		FastqWorkBuffersSE decompBuffers;
		std::vector<FastqRecord> decompReads;
#endif


#if (DEV_DEBUG_MODE)
			if (verboseMode_) std::cout << "Processing standard bins...\n" << std::flush;
#endif

		// process std bins
		//
		const uint32 totalPartsCount = extractor->GetBlockDescriptors(true).size();
		uint32 partsProcessed = 0;
		while (extractor->ExtractNextStdBin(binBin))
		{
			const uint32 minId = binBin.signature;

			ASSERT(minId != 0);
			ASSERT(binBin.metaSize != 0);

			reads.clear();
			workBuffers.Reset();
			mainPackCtx.Clear();

#if EXTRA_MEM_OPT
			reads.shrink_to_fit();
#endif
			packer.UnpackFromBin(binBin, reads, *mainPackCtx.graph,
								 mainPackCtx.stats, tmpChunk, false);

			ASSERT(binBin.rawDnaSize > 0);
			compBin.Reset();

#if (DEV_DEBUG_MODE)
			if (verboseMode_)
			{
				std::cout << minId << std::flush;
			}

			compressor.Compress(reads, mainPackCtx, minId, binBin.rawDnaSize,
								   workBuffers.fastqWorkBin, compBin);

			// decompress
			if (verboseMode_)
			{
				std::cout << "." << std::flush;
			}
			decompReads.clear();
			decompBuffers.Reset();
			decompressor.Decompress(compBin, decompReads, decompBuffers.fastqWorkBin,
									   decompBuffers.fastqBuffer);

			if (verboseMode_)
			{
				std::cout << ".\n" << std::flush;
			}

			// compare and validate
			//packer.UnpackFromBin(binBin, reads, workBuffers.fastqBuffer);
			ASSERT(reads.size() == decompReads.size());
			FastqComparator comparator;
			std::sort(reads.begin(), reads.end(), comparator);
			std::sort(decompReads.begin(), decompReads.end(), comparator);

			// validate the results
			for (uint64 i = 0; i < reads.size(); ++i)
			{
				FastqRecord& r_d = reads[i];
				FastqRecord& r_o = decompReads[i];

				//ASSERT(r_d.minimPos == r_o.minimPos);
				//ASSERT(r_d.IsReadReverse() == r_o.IsReadReverse());
				ASSERT(std::equal(r_d.seq, r_d.seq + r_d.seqLen, r_o.seq));
			}

#else
			compressor.Compress(reads, mainPackCtx, minId, binBin.rawDnaSize,
								   workBuffers.fastqWorkBin, compBin);
#endif

			stats.Update(compBin.stats);

			dnarch->WriteNextBin(compBin.dataBuffer, minId);


			if (verboseMode_)
			{
				partsProcessed++;

				std::cerr << '\r' << "Parts processed: " << partsProcessed;

				if (totalPartsCount > 0)
					std::cerr << " (" << partsProcessed * 100 / totalPartsCount << "%) ";
				std::cerr << std::flush;
			}
		}
	}

	// print stats
	//
	if (verboseMode_)
	{
		const auto streamNames = FastqWorkBuffersSE::GetBufferNames();
		//ASSERT(stats.bufferSizes.size() == streamNames.size());

		std::cout << '\n';
		if (stats.bufferSizes.count("CompSize"))
		{
			std::cout << "StreamSizes:\n";
			const auto& buf = stats.bufferSizes.at("CompSize");
			for (uint32 i = 0; i < buf.size(); ++i)
				std::cout << streamNames[i] << " " << buf[i] << '\n';
			if (stats.counts.count("NDnaCompSize"))
			{
				std::cout << "NDna: " << stats.counts.at("NDnaCompSize") << '\n';
				std::cout << "NQua: " << stats.counts.at("NQuaCompSize") << '\n';
				std::cout << "NReadIdToken: " << stats.counts.at("NReadIdTokenCompSize") << '\n';
				std::cout << "NReadIdValue: " << stats.counts.at("NReadIdValueCompSize") << '\n';
			}
			else
			{
				std::cout << "NDna: 0\n";
				std::cout << "NQua: 0\n";
				std::cout << "NReadIdToken: 0\n";
				std::cout << "NReadIdValue: 0\n";
			}

			std::cout << std::endl;
		}

		std::cout << "**** **** **** ****\n";



#if DEV_DEBUG_MODE
		std::cout << '\n';
		if (stats.bufferSizes.count("RawSize"))
		{
			std::cout << "RawStreamSizes:\n";
			const auto& buf = stats.bufferSizes.at("RawSize");
			for (uint32 i = 0; i < buf.size(); ++i)
				std::cout << streamNames[i] << " " << buf[i] << '\n';
			if (stats.counts.count("NDnaRawSize"))
			{
				std::cout << "NDna: " << stats.counts.at("NDnaRawSize") << '\n';
				std::cout << "NQua: " << stats.counts.at("NDnaRawSize") << '\n';
				std::cout << "NReadId: " << stats.counts.at("NReadIdRawSize") << '\n';
			}
			else
			{
				std::cout << "NQua: 0\n";
				std::cout << "NQua: 0\n";
				std::cout << "NReadId: 0\n";
			}

			std::cout << std::endl;
		}

		std::cout << "**** **** **** ****\n";

		const auto filterCounts = {"NDnaCompSize", "NQuaCompSize", "NReadIdTokenCompSize", "NReadIdValueCompSize", "NDnaRawSize", "NReadIdRawSize"};
		for (auto& lp : stats.counts)
		{
			if (std::find(filterCounts.begin(), filterCounts.end(), lp.first) != filterCounts.end())
				continue;

			std::cout << lp.first << ": " << lp.second << '\n';
		}

		std::cout << "**** **** **** ****\n";

		for (auto& lp : stats.freqs)
		{
			std::cout << lp.first << '\n';
			for (auto v : lp.second)
			{
				std::cout << v.first << " " << v.second << '\n';
			}

			std::cout << "**** **** **** ****\n";
		}
#endif

		std::cout << std::flush;
	}


	// store quality and headers compression data and finalize the compression
	//
	dnarch->SetQualityCompressionData(globalQuaData);

	if (binConf.archiveType.readsHaveHeaders)
	{
		dnarch->SetHeadersCompressionData(headData);
	}

	extractor->FinishDecompress();
	dnarch->FinishCompress();


	// HACK!!!
	QvzCodebook* q = (QvzCodebook *)(&extractor->GetFileFooter().quaData.codebook);
	q->qlist = NULL;

	delete dnarch;
	delete extractor;
}


void CompressorModuleSE::Dnarch2Dna(const std::string &inArchiveFile_,
								const std::string &outDnaFile_,
								uint32 threadsNum_)
{
	ArchiveFileReader* dnarch = new ArchiveFileReader();
	ArchiveFileReader::ArchiveConfig archConfig;

	dnarch->StartDecompress(inArchiveFile_, archConfig);
	ASSERT(archConfig.archType.readType == ArchiveType::READ_SE);

	CompressorParams compParams;
	compParams.archType = archConfig.archType;
	compParams.minimizer = archConfig.minParams;
	compParams.quality = archConfig.quaParams;

	FastqFileWriterSE* dnaFile = new FastqFileWriterSE(outDnaFile_);

	// here's the essential QVZ data to decompress qualities
	// (can be empty, depending on the scheme)
	//
	const QualityCompressionData& globalQuaData = dnarch->GetQualityCompressionData();
	//
	// ...

	const auto& headData = dnarch->GetHeadersCompressionData();

	if (threadsNum_ > 1)
	{
		const uint32 partNum = threadsNum_ + (threadsNum_ >> 2);
		const uint64 inBufferSize = 1 << 20;
		const uint64 outBufferSize = 1 << 20;

		typedef TRawDnaPartsWriter<FastqChunkCollectionSE> RawDnaPartsWriter;
		typedef TDnaPartsDecompressor<FastqChunkCollectionSE> DnaPartsDecompressor;
		typedef TDataPool<FastqChunkCollectionSE> FastqPartsPool;
		typedef TDataQueue<FastqChunkCollectionSE> FastqPartsQueue;

		CompressedFastqBlockPool* inPool = new CompressedFastqBlockPool(partNum, inBufferSize);
		CompressedFastqBlockQueue* inQueue = new CompressedFastqBlockQueue(partNum, 1);

		FastqPartsPool* outPool = new FastqPartsPool(partNum, outBufferSize);
		FastqPartsQueue* outQueue = new FastqPartsQueue(partNum, threadsNum_);

		ArchivePartsReader* inReader = new ArchivePartsReader(dnarch, inQueue, inPool);
		RawDnaPartsWriter* outWriter = new RawDnaPartsWriter(dnaFile, outQueue, outPool);

		// launch stuff
		//
		mt::thread readerThread(mt::ref(*inReader));

		std::vector<IOperator*> operators;

#ifdef USE_BOOST_THREAD
		boost::thread_group opThreadGroup;

		for (uint32 i = 0; i < threadsNum_; ++i)
		{
			IOperator* op = new DnaPartsDecompressor(compParams,
													 inQueue, inPool,
													 outQueue, outPool);
			operators.push_back(op);
			opThreadGroup.create_thread(mt::ref(*op));
		}

		(*outWriter)();

		readerThread.join();
		opThreadGroup.join_all();

#else
		std::vector<mt::thread> opThreadGroup;

		for (uint32 i = 0; i < threadsNum_; ++i)
		{
			IOperator* op = new DnaPartsDecompressor(compParams, globalQuaData, headData,
													 inQueue, inPool, outQueue, outPool);
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
		FastqDecompressor compressor(compParams, globalQuaData, headData);
		CompressedFastqBlock compBlock;
		FastqWorkBuffersSE  workBuffers;
		std::vector<FastqRecord> reads;
		FastqChunkCollectionSE dnaChunk;

		uint32 signatureId = 0;
		std::string signature;
		signature.resize(compParams.minimizer.signatureLen, 'N');

		while (dnarch->ReadNextBin(compBlock.dataBuffer, signatureId))
		{
			compBlock.signatureId = signatureId;
			compressor.Decompress(compBlock, reads, workBuffers.fastqWorkBin,
									 workBuffers.fastqBuffer);

			compParams.minimizer.GenerateMinimizer(signatureId, (char*)signature.c_str());
			FastqRecordsParserDynSE parser(compParams.archType.readsHaveHeaders, signature);
			parser.ParseTo(reads, dnaChunk, 1);

			dnaFile->WriteNextChunk(dnaChunk);
		}
	}

	dnarch->FinishDecompress();
	dnaFile->Close();

	delete dnaFile;
	delete dnarch;
}


void CompressorModulePE::Bin2Dnarch(const std::string &inBinFile_, const std::string &outArchiveFile_,
								const CompressorParams& compParams_, const CompressorAuxParams& auxParams_,
								uint32 threadsNum_, bool verboseMode_)
{

	BinModuleConfig binConf;
	BinFileExtractor* extractor = new BinFileExtractor(compParams_.extractor.minBinSize);

	extractor->StartDecompress(inBinFile_, binConf);

	// TODO: handle this in parent module
	//
	ASSERT(binConf.archiveType.readType == ArchiveType::READ_PE);

	// WARN: remember about synchronisation of those parameters
	//
	ArchiveFileWriter::ArchiveConfig archConfig;
	archConfig.archType = binConf.archiveType;
	archConfig.minParams = binConf.minimizer;
	archConfig.quaParams = binConf.quaParams;

	CompressorParams params = compParams_;
	params.archType = binConf.archiveType;
	params.minimizer = binConf.minimizer;
	params.quality = binConf.quaParams;
	//
	// //

	ArchiveFileWriter* dnarch = new ArchiveFileWriter();
	dnarch->StartCompress(outArchiveFile_, archConfig);

	const uint32 totalBinsCount = extractor->GetBlockDescriptors(true).size();
	uint64 nRecordsCount = 0;


	// here we can already calculate global QVZ codebooks
	//
	const QualityCompressionData& globalQuaData = extractor->GetFileFooter().quaData;

	// ...

	const auto& headData = extractor->GetFileFooter().headData;




#if (DEV_DEBUG_MODE)
	if (verboseMode_) std::cout << "Processing small bins...\n" << std::flush;
#endif

	CompressedFastqBlockStats stats;
	{
		FastqCompressor compressor(params, globalQuaData, headData, auxParams_);
		FastqNodesPackerPE packer(binConf);

#if(DEV_DEBUG_MODE)
		FastqCompressor decompressor(params, globalQuaData);
		FastqWorkBuffersPE decompBuffers;
		std::vector<FastqRecord> decompReads;
#endif

		// TODO: integrate 3 of those
		std::vector<FastqRecord> reads;
		PackContext mainPackCtx;
		FastqWorkBuffersPE workBuffers;

		CompressedFastqBlock compBin;
		BinaryBinBlock binBin;

		uint64 totalDnaBufferSize = 0;
		uint64 totalHeadBufferSize = 0;
		auto descriptors = extractor->GetBlockDescriptors(false);
		for (const auto& desc : descriptors)
		{
			totalDnaBufferSize += desc.second->totalRawDnaSize;
			totalHeadBufferSize += desc.second->totalRawHeadSize;
			nRecordsCount += desc.second->totalRecordsCount;
		}

		auto nBinDesc = extractor->GetNBlockDescriptor();
		if (nBinDesc.second != NULL)
		{
			totalDnaBufferSize += nBinDesc.second->totalRawDnaSize;
			totalHeadBufferSize += nBinDesc.second->totalRawHeadSize;
			nRecordsCount += nBinDesc.second->totalRecordsCount;
		}

		const uint64 bufferPreallocSize = totalDnaBufferSize * 2
				+ (uint64)binConf.archiveType.readsHaveHeaders * totalHeadBufferSize;
		if (workBuffers.fastqBuffer.data.Size() < bufferPreallocSize)		// support for quality
			workBuffers.fastqBuffer.data.Extend(bufferPreallocSize);

		// WARN: we need to reserve the memory for reads as
		// later we'll be using pointers to records while depacking nodes
		if (reads.size() < nRecordsCount)
			reads.reserve(nRecordsCount);


		// extract and unpack small bins
		//
		while (extractor->ExtractNextSmallBin(binBin))
		{
			ASSERT(binBin.metaSize != 0);

			packer.UnpackFromBin(binBin, reads, *mainPackCtx.graph, mainPackCtx.stats, workBuffers.fastqBuffer, true);
		}

		if (extractor->ExtractNBin(binBin))
		{
			packer.UnpackFromBin(binBin, reads, *mainPackCtx.graph, mainPackCtx.stats, workBuffers.fastqBuffer, true);
		}


		// un-reverse-compliment records
		//
		if (reads.size() > 0)
		{
			FastqRecordBuffer rcRec;

			for (FastqRecord& rec : reads)
			{
				// handle PE
				//
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


			// compress all bins together
			//
			const uint32 nSignature = params.minimizer.SignatureN();

#if (DEV_DEBUG_MODE)
			compressor.Compress(reads, mainPackCtx, nSignature,
								   totalDnaBufferSize, workBuffers.fastqWorkBin,
								   compBin);
			// checkup
			decompReads.clear();
			decompBuffers.Reset();
			decompressor.Decompress(compBin, decompReads, decompBuffers.fastqWorkBin,
									   decompBuffers.fastqBuffer);

			// compare and validate
			ASSERT(reads.size() == decompReads.size());
			FastqComparator comparator;
			std::sort(reads.begin(), reads.end(), comparator);
			std::sort(decompReads.begin(), decompReads.end(), comparator);

			// validate the results
			for (uint64 i = 0; i < reads.size(); ++i)
			{
				FastqRecord& r_d = reads[i];
				FastqRecord& r_o = decompReads[i];

				//ASSERT(r_d.minimPos == r_o.minimPos);
				//ASSERT(r_d.IsReadReverse() == r_o.IsReadReverse());
				ASSERT(std::equal(r_d.seq, r_d.seq + r_d.seqLen + r_d.auxLen, r_o.seq));
			}

#else
			compressor.Compress(reads, mainPackCtx, nSignature,
								   totalDnaBufferSize, workBuffers.fastqWorkBin,
								   compBin);
#endif

			dnarch->WriteNextBin(compBin.dataBuffer, nSignature);

			stats = compBin.stats;
		}
	}

#if 1
	if (threadsNum_ > 1)
	{
		const uint32 partNum = threadsNum_ + (threadsNum_ >> 2);
		const uint64 dnaBufferSize = 1 << 20;
		const uint64 outBufferSize = 1 << 20;


		MinimizerPartsPool* inPool = new MinimizerPartsPool(partNum, dnaBufferSize);
		MinimizerPartsQueue* inQueue = new MinimizerPartsQueue(partNum, 1);

		CompressedFastqBlockPool* outPool = new CompressedFastqBlockPool(partNum, outBufferSize);
		CompressedFastqBlockQueue* outQueue = new CompressedFastqBlockQueue(partNum, threadsNum_);

		BinPartsExtractor* inReader = new BinPartsExtractor(extractor, inQueue, inPool);

		ArchivePartsWriter* outWriter = new ArchivePartsWriter(dnarch,
															 outQueue,
															 outPool,
															 verboseMode_,
															 totalBinsCount);

		// update stats from preprocessing
		//
		outWriter->GetStats().Update(stats);


		// TODO : reuse the first part of the buffered data from pool
		//


		// launch stuff
		//
		mt::thread readerThread(mt::ref(*inReader));

		std::vector<IOperator*> operators;
		operators.resize(threadsNum_);

#ifdef USE_BOOST_THREAD
		boost::thread_group opThreadGroup;

		for (uint32 i = 0; i < threadsNum_; ++i)
		{
			operators[i] = new BinPartsCompressor(params, inQueue, inPool, outQueue, outPool);
			opThreadGroup.create_thread(mt::ref(*operators[i]));
		}

		(*outWriter)();

		readerThread.join();
		opThreadGroup.join_all();


#else
		std::vector<mt::thread> opThreadGroup;

		for (IOperator*& op : operators)
		{
			op = new BinPartsCompressor(params, auxParams_, binConf, globalQuaData, headData,
										inQueue, inPool, outQueue, outPool);
			opThreadGroup.push_back(mt::thread(mt::ref(*op)));
		}

		(*outWriter)();

		readerThread.join();

		for (mt::thread& t : opThreadGroup)
		{
			t.join();
		}

#endif
		stats = outWriter->GetStats();

		for (uint32 i = 0; i < threadsNum_; ++i)
		{
			delete operators[i];
		}

		TFREE(outWriter);
		TFREE(inReader);

		TFREE(outQueue);
		TFREE(outPool);
		TFREE(inQueue);
		TFREE(inPool);
	}
	else
#endif
	{
		FastqNodesPackerPE packer(binConf);
		FastqWorkBuffersPE workBuffers;
		FastqCompressor compressor(params, globalQuaData, headData, auxParams_);

		CompressedFastqBlock compBin;
		std::vector<FastqRecord> reads;
		PackContext mainPackCtx;

		BinaryBinBlock binBin;

#if(DEV_DEBUG_MODE)
		FastqCompressor decompressor(params, globalQuaData);
		FastqWorkBuffersPE decompBuffers;
		std::vector<FastqRecord> decompReads;
#endif


#if (DEV_DEBUG_MODE)
		if (verboseMode_) std::cout << "Processing standard bins...\n" << std::flush;
#endif

		// process std bins
		//
		const uint32 totalPartsCount = extractor->GetBlockDescriptors(true).size();
		uint32 partsProcessed = 0;
		while (extractor->ExtractNextStdBin(binBin))
		{
			const uint32 minId = binBin.signature;

			ASSERT(minId != 0);
			ASSERT(binBin.metaSize != 0);

			reads.clear();
			workBuffers.Reset();
			mainPackCtx.Clear();

#if EXTRA_MEM_OPT
			reads.shrink_to_fit();
#endif

			packer.UnpackFromBin(binBin, reads, *mainPackCtx.graph,
								 mainPackCtx.stats, workBuffers.fastqBuffer, false);

			ASSERT(binBin.rawDnaSize > 0);
			compBin.Reset();

#if (DEV_DEBUG_MODE)
			if (verboseMode_)
			{
				std::cout << minId << std::flush;
			}

			compressor.Compress(reads, mainPackCtx, minId, binBin.rawDnaSize,
								   workBuffers.fastqWorkBin, compBin);

			// decompress
			if (verboseMode_)
			{
				std::cout << "." << std::flush;
			}
			decompReads.clear();
			decompBuffers.Reset();
			decompressor.Decompress(compBin, decompReads, decompBuffers.fastqWorkBin,
									   decompBuffers.fastqBuffer);

			if (verboseMode_)
			{
				std::cout << ".\n" << std::flush;
			}

			// compare and validate
			ASSERT(reads.size() == decompReads.size());
			FastqComparator comparator;
			std::sort(reads.begin(), reads.end(), comparator);
			std::sort(decompReads.begin(), decompReads.end(), comparator);

			// validate the results
			for (uint64 i = 0; i < reads.size(); ++i)
			{
				FastqRecord& r_d = reads[i];
				FastqRecord& r_o = decompReads[i];

				//ASSERT(r_d.minimPos == r_o.minimPos);
				//ASSERT(r_d.IsReadReverse() == r_o.IsReadReverse());
				ASSERT(std::equal(r_d.seq, r_d.seq + r_d.seqLen + r_d.auxLen, r_o.seq));
			}

#else
			compressor.Compress(reads, mainPackCtx, minId, binBin.rawDnaSize,
								   workBuffers.fastqWorkBin, compBin);
#endif

			stats.Update(compBin.stats);

			dnarch->WriteNextBin(compBin.dataBuffer, minId);


			if (verboseMode_)
			{
				partsProcessed++;

				std::cerr << '\r' << "Parts processed: " << partsProcessed;

				if (totalPartsCount > 0)
					std::cerr << " (" << partsProcessed * 100 / totalPartsCount << "%) ";
				std::cerr << std::flush;
			}
		}
	}

	// print stats
	//
	if (verboseMode_)
	{
		const auto streamNames = FastqWorkBuffersPE::GetBufferNames();
		//ASSERT(stats.bufferSizes.size() == streamNames.size());

		std::cout << '\n';
		if (stats.bufferSizes.count("CompSize"))
		{
			std::cout << "StreamSizes:\n";
			const auto& buf = stats.bufferSizes.at("CompSize");
			for (uint32 i = 0; i < buf.size(); ++i)
				std::cout << streamNames[i] << " " << buf[i] << '\n';
			if (stats.counts.count("NDnaCompSize"))
			{
				std::cout << "NDna: " << stats.counts.at("NDnaCompSize") << '\n';
				std::cout << "NQua: " << stats.counts.at("NQuaCompSize") << '\n';
				std::cout << "NReadIdToken: " << stats.counts.at("NReadIdTokenCompSize") << '\n';
				std::cout << "NReadIdValue: " << stats.counts.at("NReadIdValueCompSize") << '\n';
			}
			else
			{
				std::cout << "NDna: 0\n";
				std::cout << "NQua: 0\n";
				std::cout << "NReadIdToken: 0\n";
				std::cout << "NReadIdValue: 0\n";
			}

			std::cout << std::endl;
		}

		std::cout << "**** **** **** ****\n";


#if DEV_DEBUG_MODE
		std::cout << '\n';
		if (stats.bufferSizes.count("RawSize"))
		{
			std::cout << "RawStreamSizes:\n";
			const auto& buf = stats.bufferSizes.at("RawSize");
			for (uint32 i = 0; i < buf.size(); ++i)
				std::cout << streamNames[i] << " " << buf[i] << '\n';
			if (stats.counts.count("NDnaRawSize"))
			{
				std::cout << "NDna: " << stats.counts.at("NDnaRawSize") << '\n';
				std::cout << "NQua: " << stats.counts.at("NDnaRawSize") << '\n';
				std::cout << "NReadId: " << stats.counts.at("NReadIdRawSize") << '\n';
			}
			else
			{
				std::cout << "NQua: 0\n";
				std::cout << "NQua: 0\n";
				std::cout << "NReadId: 0\n";
			}

			std::cout << std::endl;
		}

		std::cout << "**** **** **** ****\n";

		const auto filterCounts = {"NDnaCompSize", "NQuaCompSize", "NReadIdTokenCompSize", "NReadIdValueCompSize", "NDnaRawSize", "NReadIdRawSize"};
		for (auto& lp : stats.counts)
		{
			if (std::find(filterCounts.begin(), filterCounts.end(), lp.first) != filterCounts.end())
				continue;

			std::cout << lp.first << ": " << lp.second << '\n';
		}

		std::cout << "**** **** **** ****\n";


		//const auto filterFreqs = {"LzId", "LzId_PE"};
		for (auto& lp : stats.freqs)
		{
			std::cout << lp.first << '\n';
			for (auto v : lp.second)
			{
				std::cout << v.first << " " << v.second << '\n';
			}

			std::cout << "**** **** **** ****\n";
		}
#endif

		std::cout << std::flush;
	}


	// store quality and header compression data and finalize the compression
	//
	dnarch->SetQualityCompressionData(globalQuaData);

	if (binConf.archiveType.readsHaveHeaders)
	{
		dnarch->SetHeadersCompressionData(headData);
	}

	extractor->FinishDecompress();
	dnarch->FinishCompress();


	// HACK!!!
	QvzCodebook* q = (QvzCodebook *)(&extractor->GetFileFooter().quaData.codebook);
	q->qlist = NULL;


	delete dnarch;
	delete extractor;
}


void CompressorModulePE::Dnarch2Dna(const std::string &inArchiveFile_, const std::string &outDnaFile1_,
								const std::string &outDnaFile2_, uint32 threadsNum_)
{
	ArchiveFileReader* dnarch = new ArchiveFileReader();

	ArchiveFileReader::ArchiveConfig archConfig;
	dnarch->StartDecompress(inArchiveFile_, archConfig);
	ASSERT(archConfig.archType.readType == ArchiveType::READ_PE);

	CompressorParams compParams;
	compParams.archType = archConfig.archType;
	compParams.minimizer = archConfig.minParams;
	compParams.quality = archConfig.quaParams;

	FastqFileWriterPE* dnaFile = new FastqFileWriterPE(outDnaFile1_, outDnaFile2_);


	// here's the essential QVZ data to decompress qualities
	// (can be empty, depending on the scheme)
	//
	const QualityCompressionData& globalQuaData = dnarch->GetQualityCompressionData();
	//
	// ...

	const auto& headData = dnarch->GetHeadersCompressionData();


	if (threadsNum_ > 1)
	{
		const uint32 partNum = threadsNum_ + (threadsNum_ >> 2);
		const uint64 inBufferSize = 1 << 8;
		const uint64 outBufferSize = 1 << 8;

		typedef TRawDnaPartsWriter<FastqChunkCollectionPE> RawDnaPartsWriter;
		typedef TDnaPartsDecompressor<FastqChunkCollectionPE> DnaPartsDecompressor;
		typedef TDataPool<FastqChunkCollectionPE> FastqPartsPool;
		typedef TDataQueue<FastqChunkCollectionPE> FastqPartsQueue;

		CompressedFastqBlockPool* inPool = new CompressedFastqBlockPool(partNum, inBufferSize);
		CompressedFastqBlockQueue* inQueue = new CompressedFastqBlockQueue(partNum, 1);

		FastqPartsPool* outPool = new FastqPartsPool(partNum, outBufferSize);
		FastqPartsQueue* outQueue = new FastqPartsQueue(partNum, threadsNum_);

		ArchivePartsReader* inReader = new ArchivePartsReader(dnarch, inQueue, inPool);
		RawDnaPartsWriter* outWriter = new RawDnaPartsWriter(dnaFile, outQueue, outPool);

		// launch stuff
		//
		mt::thread readerThread(mt::ref(*inReader));

		std::vector<IOperator*> operators;

		std::vector<mt::thread> opThreadGroup;

		for (uint32 i = 0; i < threadsNum_; ++i)
		{
			IOperator* op = new DnaPartsDecompressor(compParams, globalQuaData, headData,
													 inQueue, inPool, outQueue, outPool);
			operators.push_back(op);
			opThreadGroup.push_back(mt::thread(mt::ref(*op)));
		}

		(*outWriter)();

		readerThread.join();

		for (mt::thread& t : opThreadGroup)
			t.join();

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
		FastqDecompressor compressor(compParams, globalQuaData, headData);
		CompressedFastqBlock compBlock;
		std::vector<FastqRecord> reads;
		FastqWorkBuffersPE  workBuffers;
		FastqChunkCollectionPE dnaChunk;

		uint32 signatureId = 0;
		std::string signature;
		signature.resize(compParams.minimizer.signatureLen, 'N');

		while (dnarch->ReadNextBin(compBlock.dataBuffer, signatureId))
		{
			compBlock.signatureId = signatureId;
			compressor.Decompress(compBlock, reads, workBuffers.fastqWorkBin,
									 workBuffers.fastqBuffer);


			compParams.minimizer.GenerateMinimizer(signatureId, (char*)signature.c_str());
			FastqRecordsParserDynPE parser(compParams.archType.readsHaveHeaders,
										headData.pairedEndFieldIdx,
										signature);
			parser.ParseTo(reads, dnaChunk, 1);

			dnaFile->WriteNextChunk(dnaChunk);
		}
	}

	dnarch->FinishDecompress();
	dnaFile->Close();

	delete dnaFile;
	delete dnarch;
}
