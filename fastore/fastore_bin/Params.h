/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_BINPARAMS
#define H_BINPARAMS

#include "Globals.h"
#include "Quality.h"

#include <algorithm>


/**
 * Parameters used across all the modules
 *
 */
struct MinimizerParameters
{
	struct Default
	{
		static const uint32 SignatureLength = 8;
		static const uint32 SkipZoneLength = 0;
		static const uint32 SignatureMaskCutoffBits = 0;
	};

	uint8 signatureLen;
	uint8 skipZoneLen;					// TODO: add to filtering opts
	uint8 signatureMaskCutoffBits;		// ---
	char dnaSymbolOrder[5];

	MinimizerParameters()
		:	signatureLen(Default::SignatureLength)
		,	skipZoneLen(Default::SignatureLength)
		,	signatureMaskCutoffBits(Default::SignatureMaskCutoffBits)
	{
		dnaSymbolOrder[0] = 'A';
		dnaSymbolOrder[1] = 'C';
		dnaSymbolOrder[2] = 'G';
		dnaSymbolOrder[3] = 'T';
		dnaSymbolOrder[4] = 'N';
	}

	MinimizerParameters(uint32 minLen_, uint32 cutoffLen_, const char* symOrder_ = "ACGTN")
		:	signatureLen(minLen_)
		,	skipZoneLen(cutoffLen_)
	{
		std::copy(symOrder_, symOrder_ + 5, dnaSymbolOrder);
	}

	uint64 TotalMinimizersCount() const
	{
		return 1 << (signatureLen * 2);
	}

	uint64 SignatureN() const
	{
		return TotalMinimizersCount();
	}

	void GenerateMinimizer(uint32 minimizerId_, char* buf_)	const // WARNING: no boundary check !!!
	{
		ASSERT(minimizerId_ <= TotalMinimizersCount());

		if (minimizerId_ == TotalMinimizersCount())
		{
			for (int32 i = 0; i < signatureLen; i++)
				buf_[i] = 'N';
		}
		else
		{
			for (int32 i = signatureLen - 1; i >= 0; --i)
			{
				buf_[i] = dnaSymbolOrder[minimizerId_ & 0x03];
				minimizerId_ >>= 2;
			}
		}
	}

	uint32 ReverseSignature(uint32 signature_) const
	{
		uint32 sigRevLookup[4] = {3, 2, 1, 0};
		uint32 revSig = 0;
		for (uint32 i = 0; i < signatureLen; ++i)
		{
			revSig <<= 2;
			revSig |= sigRevLookup[signature_ & 3];
			signature_ >>= 2;
		}
		return revSig;
	}
};


struct MinimizerFilteringParameters
{
	struct Default
	{
		static const bool FilterLowQualitySignatures = false;
		static const uint8 LowQualityThresholdValue = 6;
	};

	bool filterLowQualitySignatures;
	uint8 lowQualityThreshold;

	MinimizerFilteringParameters()
		:	filterLowQualitySignatures(Default::FilterLowQualitySignatures)
		,	lowQualityThreshold(Default::LowQualityThresholdValue)
	{}
};


struct CategorizerParameters
{
	static const uint32 DefaultMinimumPartialBinSize = 8;

	uint32 minBlockBinSize;

	CategorizerParameters()
		:	minBlockBinSize(DefaultMinimumPartialBinSize)
	{}
};


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
