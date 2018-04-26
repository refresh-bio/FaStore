/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_BINBLOCKDATA
#define H_BINBLOCKDATA

#include "Globals.h"

#include <vector>
#include <map>

#include "../core/Buffer.h"
#include "../qvz/Stats.h"


/**
 * Descriptor for a partial (or full) bin
 * containing DNA, ID and QUA data
 *
 */
struct BinaryBinDescriptor
{
	uint64 metaSize;
	uint64 dnaSize;
	uint64 quaSize;
	uint64 headSize;

	uint64 recordsCount;
	uint64 rawDnaSize;
	uint64 rawHeadSize;

	BinaryBinDescriptor()
	{
		Clear();
	}

	BinaryBinDescriptor(const BinaryBinDescriptor& b_) = default;

	void Clear()
	{
		metaSize = 0;
		dnaSize = 0;
		quaSize = 0;
		headSize = 0;
		recordsCount = 0;
		rawDnaSize = 0;
		rawHeadSize = 0;
	}
};


/**
 * Stores all the bin descriptors to read/write data from disk
 *
 */
struct BinaryBinBlock
{
	enum BlockTypeEnum
	{
		MultiSignatureType,
		SingleSignatureType
	};

	static const uint64 DefaultMetaBufferSize = 1 << 6;
	static const uint64 DefaultDnaBufferSize = 1 << 8;
	static const uint64 DefaultQualityBufferSize = 1 << 8;
	static const uint64 DefaultHeaderBufferSize = 16;

	// TODO: refactor
	BlockTypeEnum blockType;

	std::map<uint32, BinaryBinDescriptor> descriptors;
	uint32 signature;

	std::vector<BinaryBinDescriptor> auxDescriptors;

	FastqRawBlockStats stats;
	//


	// TODO: vector of buffers
	//
	Buffer metaData;
	Buffer dnaData;
	Buffer quaData;
	Buffer headData;		// optional

	uint64 metaSize;
	uint64 dnaSize;
	uint64 quaSize;
	uint64 headSize;		// optional
	//
	// //

	uint64 rawDnaSize;		// == raw quality size
	uint64 rawHeadSize;		// optional

	BinaryBinBlock(uint64 dnaBufferSize_ = DefaultDnaBufferSize,
				   uint64 metaBufferSize_ = DefaultMetaBufferSize,
				   uint64 quaBufferSize_ = DefaultQualityBufferSize,
				   uint64 headBufferSize_ = DefaultHeaderBufferSize)
		:	blockType(MultiSignatureType)
		,	signature(0)
		,	metaData(metaBufferSize_)
		,	dnaData(dnaBufferSize_)
		,	quaData(quaBufferSize_)
		,	headData(headBufferSize_)
	{
		Clear();
	}

	void Clear()
	{
		signature = 0;

		metaSize = 0;
		dnaSize = 0;
		quaSize = 0;

		rawDnaSize = 0;

		headSize = 0;
		rawHeadSize = 0;

		descriptors.clear();
		auxDescriptors.clear();

		stats.Clear();

#if EXTRA_MEM_OPT
		if (auxDescriptors.capacity() > 0)
			auxDescriptors.shrink_to_fit();

		if (metaData.Size() > DefaultMetaBufferSize)
			metaData.Shrink(DefaultMetaBufferSize);

		if (dnaData.Size() > DefaultDnaBufferSize)
			dnaData.Shrink(DefaultDnaBufferSize);

		if (quaData.Size() > DefaultQualityBufferSize)
			quaData.Shrink(DefaultQualityBufferSize);

		if (headData.Size() > DefaultHeaderBufferSize)
			headData.Shrink(DefaultHeaderBufferSize);
#endif
	}

	void Reset()		// for compatibility with queues, TODO: unify with 'Clear()'
	{
		Clear();
	}

	void Swap(BinaryBinBlock& b_)
	{
		std::swap(blockType, b_.blockType);
		std::swap(signature, b_.signature);
		std::swap(metaSize, b_.metaSize);
		std::swap(dnaSize, b_.dnaSize);
		std::swap(quaSize, b_.quaSize);
		std::swap(headSize, b_.headSize);
		std::swap(rawDnaSize, b_.rawDnaSize);
		std::swap(rawHeadSize, b_.rawHeadSize);

		descriptors.swap(b_.descriptors);
		auxDescriptors.swap(b_.auxDescriptors);

		metaData.Swap(b_.metaData);
		dnaData.Swap(b_.dnaData);
		quaData.Swap(b_.quaData);
		headData.Swap(b_.headData);

		std::swap(stats, b_.stats);
	}
};


#endif // H_BINBLOCKDATA
