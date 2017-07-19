/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_DNARCHOPERATOR
#define H_DNARCHOPERATOR

#include "../fastore_bin/Globals.h"
#include "../fastore_bin/DataPool.h"
#include "../fastore_bin/DataQueue.h"
#include "../fastore_bin/FastqPacker.h"
#include "../fastore_bin/FastqParser.h"

#include "../fastore_rebin/Params.h"
#include "../fastore_rebin/RebinOperator.h"

#include "Params.h"
#include "CompressedBlockData.h"
#include "BinFileExtractor.h"
#include "ArchiveFile.h"
#include "FastqCompressor.h"
#include "CompressorOperator.h"

#include "../fastore_bin/QVZ.h"


typedef TDataQueue<CompressedFastqBlock> CompressedFastqBlockQueue;
typedef TDataPool<CompressedFastqBlock> CompressedFastqBlockPool;



/**
 * Compresses binned reads block-wise.
 * Used in multi-threaded mode.
 *
 */
class BinPartsCompressor : public IOperator
{
public:
	BinPartsCompressor(const CompressorParams& compParams_, const CompressorAuxParams& auxCompParams_,
					   const BinModuleConfig& binConf_, const QualityCompressionData& globalQuaData_,
					   const FastqRawBlockStats::HeaderStats& headerData_,
						MinimizerPartsQueue* inPartsQueue_, MinimizerPartsPool* inPartsPool_,
						CompressedFastqBlockQueue* outPartsQueue_, CompressedFastqBlockPool* outPartsPool_)
		:	compParams(compParams_)
		,	headerData(headerData_)
		,	binConf(binConf_)
		,	globalQuaData(globalQuaData_)
		,	auxCompParams(auxCompParams_)
		,	inPartsQueue(inPartsQueue_)
		,	inPartsPool(inPartsPool_)
		,	outPartsQueue(outPartsQueue_)
		,	outPartsPool(outPartsPool_)
	{}

	void Run();

private:
	const CompressorParams compParams;
	const FastqRawBlockStats::HeaderStats& headerData;

	const BinModuleConfig& binConf;
	const QualityCompressionData& globalQuaData;
	const CompressorAuxParams auxCompParams;

	MinimizerPartsQueue* inPartsQueue;
	MinimizerPartsPool* inPartsPool;
	CompressedFastqBlockQueue* outPartsQueue;
	CompressedFastqBlockPool* outPartsPool;
};



/**
 * Writes to a given stream compressed FASTQ reads in blocks popped from the queue.
 * Used in multi-threaded mode.
 *
 */
class ArchivePartsWriter : public IOperator
{
public:
	ArchivePartsWriter(ArchiveFileWriter* partsStream_,
					   CompressedFastqBlockQueue* partsQueue_,
					   CompressedFastqBlockPool* partsPool_,
					   bool verboseMode_ = false,
					   uint32 totalPartsCount_ = 0)
		:	verboseMode(verboseMode_)
		,	totalPartsCount(totalPartsCount_)
		,	partsStream(partsStream_)
		,	partsQueue(partsQueue_)
		,	partsPool(partsPool_)
	{}

	void Run();

	const CompressedFastqBlockStats& GetStats() const
    {
        return stats;
    }

	CompressedFastqBlockStats& GetStats()
    {
        return stats;
    }

private:
	const bool verboseMode;
	const uint32 totalPartsCount;

	ArchiveFileWriter* partsStream;
	CompressedFastqBlockQueue* partsQueue;
	CompressedFastqBlockPool* partsPool;

	CompressedFastqBlockStats stats;
};


/**
 * Reads from a given stream compressed FASTQ reads in blocks pushing them to a queue.
 * Used in multi-threaded mode.
 *
 */
class ArchivePartsReader : public IOperator
{
public:
	ArchivePartsReader(ArchiveFileReader* partsStream_,
					   CompressedFastqBlockQueue* partsQueue_,
					   CompressedFastqBlockPool* partsPool_)
		:	partsStream(partsStream_)
		,	partsQueue(partsQueue_)
		,	partsPool(partsPool_)
	{}

	void Run();

private:
	ArchiveFileReader* partsStream;
	CompressedFastqBlockQueue* partsQueue;
	CompressedFastqBlockPool* partsPool;
};


/**
 * Decompresses FASTQ reads in blocks and pushes them to a queue.
 * Used in multi-threaded mode.
 *
 */
template <class _TChunkType>
class TDnaPartsDecompressor : public IOperator
{
public:
	typedef _TChunkType ChunkType;
	typedef TDataPool<ChunkType> FastqPartsPool;
	typedef TDataQueue<ChunkType> FastqPartsQueue;

	TDnaPartsDecompressor(const CompressorParams& compParams_,
						  const QualityCompressionData& globalQuaData_,
						  const FastqRawBlockStats::HeaderStats& headerData_,
						  CompressedFastqBlockQueue* inPartsQueue_,
						  CompressedFastqBlockPool* inPartsPool_,
						  FastqPartsQueue* outPartsQueue_, FastqPartsPool* outPartsPool_,
						  uint64 readIdxOffset_ = 0)
		:	compParams(compParams_)
		,	globalQuaData(globalQuaData_)
		,	headerData(headerData_)
		,	inPartsQueue(inPartsQueue_)
		,	inPartsPool(inPartsPool_)
		,	outPartsQueue(outPartsQueue_)
		,	outPartsPool(outPartsPool_)
	{
		(void)readIdxOffset_;
	}

	void Run();

protected:
	const CompressorParams compParams;
	const QualityCompressionData& globalQuaData;
	const FastqRawBlockStats::HeaderStats& headerData;

	CompressedFastqBlockQueue* inPartsQueue;
	CompressedFastqBlockPool* inPartsPool;
	FastqPartsQueue* outPartsQueue;
	FastqPartsPool* outPartsPool;
};


/**
 * Writes to stream decompressed FASTQ reads in blocks popped from a queue.
 * Used in multi-threaded mode.
 *
 */
template <class _TChunkType>
class TRawDnaPartsWriter : public IOperator
{
public:
	typedef _TChunkType ChunkType;
	typedef TDataPool<ChunkType> FastqPartsPool;
	typedef TDataQueue<ChunkType> FastqPartsQueue;

	TRawDnaPartsWriter(IFastqStreamWriter* partsStream_, FastqPartsQueue* partsQueue_, FastqPartsPool* partsPool_)
		:	partsStream(partsStream_)
		,	partsQueue(partsQueue_)
		,	partsPool(partsPool_)
	{}

	void Run();

private:
	IFastqStreamWriter* partsStream;
	FastqPartsQueue* partsQueue;
	FastqPartsPool* partsPool;
};


#endif // H_DNARCHOPERATOR
