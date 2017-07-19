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
#include "Params.h"


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
