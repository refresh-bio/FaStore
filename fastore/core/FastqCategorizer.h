/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_DNACATEGORIZER
#define H_DNACATEGORIZER

#include "Globals.h"

#include <vector>
#include <map>

#include "FastqRecord.h"


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


/**
 * Distributes FASTQ reads into bins according to their signatures
 *
 */
class FastqCategorizerBase
{
public:
	FastqCategorizerBase(const MinimizerParameters& params_,
						 const MinimizerFilteringParameters& filter_ = MinimizerFilteringParameters(),
						 const CategorizerParameters& catParams_ = CategorizerParameters());

	virtual ~FastqCategorizerBase() {}

	std::pair<uint32, uint16> FindMinimizer(const FastqRecord& rec_);
	std::map<uint32, uint16> FindMinimizers(const FastqRecord &rec_, uint32 startOff_ = 0, uint32 endCutoff_ = 0);

protected:
	const MinimizerParameters& params;
	const MinimizerFilteringParameters& filter;
	const CategorizerParameters catParams;

	const uint32 maxLongMinimValue;
	const uint32 nBinValue;

	std::array<char, 128> symbolIdxTable;
	std::vector<bool> validBinSignatures;

	uint32 ComputeMinimizer(const char* dna_, uint32 mLen_);
	bool IsMinimizerValid(uint32 minim_, uint32 mLen_);
	bool IsMinimizerQualityValid(const char* qua_, uint32 mLen_, uint32 minQ_);

	void InitializeValidBinSignatures(uint32 signatureLen_);
};


// TODO: templatize and reuse the code
//
class FastqCategorizerSE : public FastqCategorizerBase
{
public:
	FastqCategorizerSE(const MinimizerParameters& params_,
					   const MinimizerFilteringParameters& filter_ = MinimizerFilteringParameters(),
					   const CategorizerParameters& catParams_ = CategorizerParameters())
		:	FastqCategorizerBase(params_, filter_, catParams_)
	{}

	void Categorize(std::vector<FastqRecord>& records_,
					std::map<uint32, FastqRecordsPtrBin>& bins_);

protected:
	virtual void DistributeToBins(std::vector<FastqRecord>& records_,
								  std::map<uint32, FastqRecordsPtrBin>& bins_);

	using FastqCategorizerBase::FindMinimizer;
	using FastqCategorizerBase::FindMinimizers;
};


class FastqCategorizerPE : public FastqCategorizerSE
{
public:
	FastqCategorizerPE(const MinimizerParameters& params_,
					   const MinimizerFilteringParameters& filter_ = MinimizerFilteringParameters(),
					   const CategorizerParameters& catParams_ = CategorizerParameters())
		:	FastqCategorizerSE(params_, filter_, catParams_)
	{}

protected:
	virtual void DistributeToBins(std::vector<FastqRecord>& records_,
								  std::map<uint32, FastqRecordsPtrBin>& bins_);
};



#endif // H_DNACATEGORIZER
