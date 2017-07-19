/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_REBINOPERATOR
#define H_REBINOPERATOR

#include "Params.h"
#include "../fastore_bin/FastqRecord.h"
#include "../fastore_bin/Globals.h"
#include "../fastore_bin/BinOperator.h"
//#include "../fastore_pack/CompressorOperator.h"
#include "DnaRebalancer.h"

#include "../fastore_pack/BinFileExtractor.h"


typedef TDataPool<BinaryBinBlock> MinimizerPartsPool;
typedef TDataQueue<BinaryBinBlock> MinimizerPartsQueue;


/**
 * Keeps the temporary buffers to perform re-binning
 */
struct RebinWorkBuffer
{
	std::vector<FastqRecord> reads;
	RebinContext rebinCtx;
	FastqChunk dataBuffer;
	std::map<uint32, MatchNodesPtrBin> nodesMap;

	void Reset()
	{
		dataBuffer.Reset();

		// TODO: this was freshly added, check me!
		reads.clear();
		rebinCtx.Clear();
		nodesMap.clear();

#if EXTRA_MEM_OPT
		reads.shrink_to_fit();

		if (dataBuffer.data.Size() > FastqChunk::DefaultBufferSize)
			dataBuffer.data.Shrink(FastqChunk::DefaultBufferSize);
#endif
	}
};


/**
 * Re-bins the reads in multi-threaded mode
 */
class BinBalancer : public IOperator
{
public:
	BinBalancer(const BinModuleConfig& binConfig_, const BinBalanceParameters balanceParams_,
				MinimizerPartsQueue* inPartsQueue_, MinimizerPartsPool* inPartsPool_,
				BinaryPartsQueue* outPartsQueue_, BinaryPartsPool* outPartsPool_)
		:	binConfig(binConfig_)
		,	balanceParams(balanceParams_)
		,	inPartsQueue(inPartsQueue_)
		,	inPartsPool(inPartsPool_)
		,	outPartsQueue(outPartsQueue_)
		,	outPartsPool(outPartsPool_)
	{}

	void Run();

private:
	const BinModuleConfig& binConfig;
	const BinBalanceParameters balanceParams;

	MinimizerPartsQueue* inPartsQueue;
	MinimizerPartsPool* inPartsPool;
	BinaryPartsQueue* outPartsQueue;
	BinaryPartsPool* outPartsPool;
};


/**
 * Extracts the bins
 */
class BinPartsExtractor : public IOperator
{
public:
	BinPartsExtractor(BinFileExtractor* partsStream_, MinimizerPartsQueue* partsQueue_, MinimizerPartsPool* partsPool_)
		:	partsStream(partsStream_)
		,	partsQueue(partsQueue_)
		,	partsPool(partsPool_)
	{}

	void Run()
	{
		int64 partId = 0;

		BinaryBinBlock* part = NULL;
		partsPool->Acquire(part);

		while (partsStream->ExtractNextStdBin(*part))
		{
			if (part->metaSize == 0)
				continue;

			ASSERT(part->signature != 0);
			partsQueue->Push(partId++, part);

			part = NULL;
			partsPool->Acquire(part);
		}
		partsPool->Release(part);

		partsQueue->SetCompleted();
	}

private:
	BinFileExtractor* partsStream;
	MinimizerPartsQueue* partsQueue;
	MinimizerPartsPool* partsPool;
};


#endif // H_BINOPERATOR
