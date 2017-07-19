/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "../fastore_bin/Globals.h"

#include "CompressorOperator.h"
#include "BinFileExtractor.h"
#include "ArchiveFile.h"
#include "FastqCompressor.h"
#include "../fastore_bin/FastqParser.h"
#include "../fastore_rebin/NodesPacker.h"

#include <iostream>
#include <memory>

void BinPartsCompressor::Run()
{
	const bool pairedEnd = compParams.archType.readType == ArchiveType::READ_PE;

	int64 partId = 0;

	std::unique_ptr<IFastqNodesPackerDyn> packer(!pairedEnd
		? (IFastqNodesPackerDyn*)(new FastqNodesPackerDynSE(binConf))
		: (IFastqNodesPackerDyn*)(new FastqNodesPackerDynPE(binConf)));
	std::unique_ptr<IFastqWorkBuffer> workBuffers(!pairedEnd
													? new FastqWorkBuffersSE()
													: new FastqWorkBuffersPE());
	FastqCompressor compressor(compParams, globalQuaData, headerData, auxCompParams);

	PackContext mainPackCtx;
	std::vector<FastqRecord> reads;

#if (DEV_DEBUG_MODE)
	FastqCompressor decompressor(compParams, globalQuaData);
	std::unique_ptr<IFastqWorkBuffer> decompBuffers(!pairedEnd
													  ? new FastqWorkBuffersSE()
													  : new FastqWorkBuffersPE());

	std::vector<FastqRecord> decompReads;
#endif

	IFastqChunkCollection tmpChunks;

	BinaryBinBlock* inPart = NULL;
	while (inPartsQueue->Pop(partId, inPart))
	{
		ASSERT(inPart->metaSize > 0);
		ASSERT(inPart->dnaSize > 0);
		ASSERT(inPart->rawDnaSize > 0);

		const uint32 signature = inPart->signature;
		ASSERT(signature != 0);

		CompressedFastqBlock* outPart = NULL;
		outPartsPool->Acquire(outPart);

		packer->UnpackFromBin(*inPart, reads, *mainPackCtx.graph,
							  mainPackCtx.stats, tmpChunks,
							  false);

		const uint64 rawDnaSize = inPart->rawDnaSize;



#if (DEV_DEBUG_MODE)

		// compress reads
		//
		compressor.Compress(reads, mainPackCtx, signature,
							   rawDnaSize, workBuffers->fastqWorkBin,
							   *outPart);

		// decompress reads
		//
		decompReads.clear();
		decompBuffers->Reset();
		decompressor.Decompress(*outPart,decompReads,
								   decompBuffers->fastqWorkBin,
								   decompBuffers->fastqBuffer);

		outPartsQueue->Push(partId, outPart);
		outPart = NULL;

		// unpack the input reads
		//
		//reads.clear();
		//workBuffers.Reset();
		//packer.UnpackFromBin(*inPart, reads, workBuffers.fastqBuffer);

		inPartsPool->Release(inPart);
		inPart = NULL;


		// compare the decompressed reads with the input ones
		//
		ASSERT(reads.size() == decompReads.size());
		FastqComparator comparator;
		std::sort(reads.begin(), reads.end(), comparator);
		std::sort(decompReads.begin(), decompReads.end(), comparator);



		// validate the reads one-by-one
		//
		for (uint64 i = 0; i < reads.size(); ++i)
		{
			FastqRecord& r_d = reads[i];
			FastqRecord& r_o = decompReads[i];

			ASSERT(std::equal(r_d.seq, r_d.seq + r_d.seqLen + r_d.auxLen, r_o.seq));
		}

		// clear buffers
		//
		workBuffers->Reset();
		mainPackCtx.Clear();
		reads.clear();

#if EXTRA_MEM_OPT
		reads.shrink_to_fit();
		tmpChunks.Clear();
#endif

#else

		// reclaim used memory from the input part
		//
		inPart->Reset();

		inPartsPool->Release(inPart);
		inPart = NULL;

		compressor.Compress(reads,
							   mainPackCtx,
							   signature,
							   rawDnaSize,
							   workBuffers->fastqWorkBin,
							   *outPart);

		// clear buffers
		//
		workBuffers->Reset();
		mainPackCtx.Clear();
		reads.clear();

#if EXTRA_MEM_OPT
		reads.shrink_to_fit();
		tmpChunks.Clear();
#endif


		outPartsQueue->Push(partId, outPart);
		outPart = NULL;

#endif
	}

	outPartsQueue->SetCompleted();
}


void ArchivePartsWriter::Run()
{
	int64 partId = 0;
	CompressedFastqBlock* part = NULL;

	uint32 partsProcessed = 0;
	while (partsQueue->Pop(partId, part))
	{
		partsStream->WriteNextBin(part->dataBuffer, part->signatureId);

        stats.Update(part->stats);

		// reclaim used memory from the part
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


void ArchivePartsReader::Run()
{
	int64 partId = 0;
	CompressedFastqBlock* part = NULL;
	partsPool->Acquire(part);

	uint32 signatureId = 0;
	while (partsStream->ReadNextBin(part->dataBuffer, signatureId))
	{
		part->signatureId = signatureId;
		ASSERT(part->dataBuffer.size > 0);		// hack

		partsQueue->Push(partId, part);
		part = NULL;

		partsPool->Acquire(part);
		partId++;
	}
	partsPool->Release(part);
	partsQueue->SetCompleted();
}


template <class _TChunkType>
void TDnaPartsDecompressor<_TChunkType>::Run()
{
	const bool pairedEnd = compParams.archType.readType == ArchiveType::READ_PE;

	FastqDecompressor compressor(compParams, globalQuaData, headerData);

	int64 partId = 0;
	CompressedFastqBlock* inPart = NULL;

	std::unique_ptr<IFastqWorkBuffer> workBuffers(!pairedEnd
												  ? new FastqWorkBuffersSE()
												  : new FastqWorkBuffersPE);
	std::vector<FastqRecord> reads;

	std::string signature;
	signature.resize(compParams.minimizer.signatureLen, 'N');

	while (inPartsQueue->Pop(partId, inPart))
	{
		ChunkType* outPart = NULL;
		outPartsPool->Acquire(outPart);

		compressor.Decompress(*inPart, reads,
								 workBuffers->fastqWorkBin,
								 workBuffers->fastqBuffer);

		compParams.minimizer.GenerateMinimizer(inPart->signatureId, (char*)signature.c_str());

		// TODO: refactor this
		//
		std::unique_ptr<IRecordsParser> parser(!pairedEnd
											 ? (IRecordsParser*)(new FastqRecordsParserDynSE(compParams.archType.readsHaveHeaders,
																						  signature))
											 : (IRecordsParser*)(new FastqRecordsParserDynPE(compParams.archType.readsHaveHeaders,
																						  headerData.pairedEndFieldIdx,
																						  signature)));


		parser->ParseTo(reads, *outPart, 1);		// TODO: refactor, skip this step

		outPartsQueue->Push(partId, outPart);


		// reclaim used space from the part
		//
		inPart->Reset();

		inPartsPool->Release(inPart);
		inPart = NULL;


		// reclaim used space from the work buffers
		//
		workBuffers->Reset();
		reads.clear();

#if EXTRA_MEM_OPT
		reads.shrink_to_fit();
#endif

	}
	outPartsQueue->SetCompleted();
}


template <class _TChunkType>
void TRawDnaPartsWriter<_TChunkType>::Run()
{
	int64 partId = 0;
	ChunkType* part = NULL;

	while (partsQueue->Pop(partId, part))
	{
		partsStream->WriteNextChunk(*part);

#if EXTRA_MEM_OPT
		part->Clear();
#endif
		partsPool->Release(part);
	}
}



template class TRawDnaPartsWriter<FastqChunkCollectionSE>;
template class TRawDnaPartsWriter<FastqChunkCollectionPE>;

template class TDnaPartsDecompressor<FastqChunkCollectionSE>;
template class TDnaPartsDecompressor<FastqChunkCollectionPE>;
