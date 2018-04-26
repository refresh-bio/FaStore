/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_BINPARAMS
#define H_BINPARAMS

#include "../core/Globals.h"
#include "../core/FastqCategorizer.h"
#include "../qvz/Quality.h"


struct ArchiveType
{
	static const uint32 StandardQualityOffset = 33;
	static const uint32 Illumina64QualityOffset = 64;

	enum ReadType
	{
		READ_SE = 0,
		READ_PE
	};

	byte readType;
	byte qualityOffset;
	bool readsHaveHeaders;


	ArchiveType()
		:	readType(READ_SE)
		,	qualityOffset(StandardQualityOffset)
		,	readsHaveHeaders(false)
	{}
};

struct HeadersCompressionParams
{
	struct Default
	{
		static const bool PreserveComments = true;
	};

	HeadersCompressionParams()
		:	preserveComments(Default::PreserveComments)
	{}

	bool preserveComments;
	//std::array[32] libraryName;
};

struct BinModuleConfig
{
	enum BinningType
	{
		BIN_RECORDS = 0,
		BIN_NODES
	};

	static const uint64 DefaultFastqBlockSize = 1 << 28;	// 256 MB

	ArchiveType archiveType;
	CategorizerParameters catParams;
	MinimizerParameters minimizer;
	MinimizerFilteringParameters minFilter;
	QualityCompressionParams quaParams;
	HeadersCompressionParams headParams;

	uint64 fastqBlockSize;
	uint32 binningLevel;
	byte binningType;

	BinModuleConfig()
		:	fastqBlockSize(DefaultFastqBlockSize)
		,	binningLevel(0)
		,	binningType(BIN_RECORDS)
	{}
};


#endif // H_BINPARAMS
