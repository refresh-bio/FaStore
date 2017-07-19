/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_BINOPERATOR
#define H_BINOPERATOR

#include "Globals.h"
#include "DataPool.h"
#include "DataQueue.h"
#include "BinBlockData.h"
#include "Params.h"
#include "FastqRecord.h"

#include <vector>
#include <map>

#include "FastqStream.h"
#include "BinFile.h"
#include "FastqParser.h"
#include "FastqCategorizer.h"
#include "FastqPacker.h"


// operators for multi threaded processing
//
typedef TDataPool<BinaryBinBlock> BinaryPartsPool;
typedef TDataQueue<BinaryBinBlock> BinaryPartsQueue;


/**
 * Reads FASTQ files chunk-wise pushing the chunks to a queue.
 * Used in multithreaded processing.
 *
 */
template <class _TStreamReader, class _TInputFastqChunk>
class TFastqChunkReader : public IOperator
{
public:
	typedef _TInputFastqChunk InFastqChunk;
	typedef TDataPool<InFastqChunk> FastqChunkPool;
	typedef TDataQueue<InFastqChunk> FastqChunkQueue;
	typedef _TStreamReader IFastqStreamReader;

	TFastqChunkReader(IFastqStreamReader* partsStream_, FastqChunkQueue* partsQueue_, FastqChunkPool* partsPool_)
		:	partsStream(partsStream_)
		,	partsQueue(partsQueue_)
		,	partsPool(partsPool_)
	{}

	void Run()
	{
		uint64 partId = 0;

		InFastqChunk* part = NULL;
		partsPool->Acquire(part);

		while (partsStream->ReadNextChunk(*part))
		{
			partsQueue->Push(partId++, part);

			partsPool->Acquire(part);
		}

		partsPool->Release(part);

		partsQueue->SetCompleted();
	}

private:
	IFastqStreamReader* partsStream;
	FastqChunkQueue* partsQueue;
	FastqChunkPool* partsPool;
};



/**
 * Writes binned FASTQ reads pushing the blocks to a queue.
 * Used in multithreaded processing.
 *
 */
class BinChunkWriter : public IOperator
{
public:
	BinChunkWriter(BinFileWriter* partsStream_,
				   BinaryPartsQueue* partsQueue_,
				   BinaryPartsPool* partsPool_,
				   bool verboseMode_ = false,
				   uint64 totalPartsCount_ = 0)
		:	verboseMode(verboseMode_)
		,	totalPartsCount(totalPartsCount_)
		,	partsStream(partsStream_)
		,	partsQueue(partsQueue_)
		,	partsPool(partsPool_)
	{}

	void Run();

private:
	const bool verboseMode;
	const uint64 totalPartsCount;

	BinFileWriter* partsStream;
	BinaryPartsQueue* partsQueue;
	BinaryPartsPool* partsPool;
};



/**
 * Bins the FASTQ reads into bins in single-end mode.
 * Used in multithreaded processing.
 *
 */
class BinEncoderSE : public IOperator
{
public:
	typedef TDataPool<FastqChunkCollectionSE> FastqChunkPool;
	typedef TDataQueue<FastqChunkCollectionSE> FastqChunkQueue;

	BinEncoderSE(const BinModuleConfig& binConfig_,
				 FastqChunkQueue* fqPartsQueue_, FastqChunkPool* fqPartsPool_,
				 BinaryPartsQueue* binPartsQueue_, BinaryPartsPool* binPartsPool_)
		:	binConfig(binConfig_)
		,	fqPartsQueue(fqPartsQueue_)
		,	fqPartsPool(fqPartsPool_)
		,	binPartsQueue(binPartsQueue_)
		,	binPartsPool(binPartsPool_)
	{}

	void Run();

private:
	const BinModuleConfig& binConfig;

	FastqChunkQueue* fqPartsQueue;
	FastqChunkPool* fqPartsPool;
	BinaryPartsQueue* binPartsQueue;
	BinaryPartsPool* binPartsPool;
};


/**
 * Bins FASTQ reads into bins in paired-end mode.
 * Used in multithreaded processing.
 *
 */
class BinEncoderPE : public IOperator
{
public:
	typedef TDataPool<FastqChunkCollectionPE> FastqChunkPool;
	typedef TDataQueue<FastqChunkCollectionPE> FastqChunkQueue;

	BinEncoderPE(const BinModuleConfig& binConfig_,
				 FastqChunkQueue* fqPartsQueue_, FastqChunkPool* fqPartsPool_,
				 BinaryPartsQueue* binPartsQueue_, BinaryPartsPool* binPartsPool_)
		:	binConfig(binConfig_)
		,	fqPartsQueue(fqPartsQueue_)
		,	fqPartsPool(fqPartsPool_)
		,	binPartsQueue(binPartsQueue_)
		,	binPartsPool(binPartsPool_)
	{}

	void Run();

private:
	const BinModuleConfig& binConfig;

	FastqChunkQueue* fqPartsQueue;
	FastqChunkPool* fqPartsPool;
	BinaryPartsQueue* binPartsQueue;
	BinaryPartsPool* binPartsPool;
};


#endif // H_BINOPERATOR
