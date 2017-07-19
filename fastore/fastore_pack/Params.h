/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_PACKPARAMS
#define H_PACKPARAMS

#include "../fastore_bin/Globals.h"
#include "../fastore_bin/Params.h"


struct ReadsClassifierParams
{
	struct Default
	{
		static const int32 MaxCostValue = (uint16)-1;
		static const int32 AutoEncodeThresholdValue = 0;
		static const int32 ShiftCost = 1;
		static const int32 MismatchCost = 2;
		static const uint32 MaxLzWindowSize = MAX_LZ_SE;
		static const uint32 MaxPairLzWindowSize = MAX_LZ_PE;
		static const bool ExtraReduceHardReads = false;			// temporary, for backwards-compatibility
		static const bool ExtraReduceExpensiveLzMatches = false;
	};

	int32 maxCostValue;
	int32 encodeThresholdValue;
	int32 pairEncodeThresholdValue;
	int32 shiftCost;
	int32 mismatchCost;
	uint32 maxLzWindowSize;
	uint32 maxPairLzWindowSize;
	bool extraReduceHardReads;
	bool extraReduceExpensiveLzMatches;

	ReadsClassifierParams()
		:	maxCostValue(Default::MaxCostValue)
		,	encodeThresholdValue(Default::AutoEncodeThresholdValue)
		,	pairEncodeThresholdValue(Default::AutoEncodeThresholdValue)
		,	shiftCost(Default::ShiftCost)
		,	mismatchCost(Default::MismatchCost)
		,	maxLzWindowSize(Default::MaxLzWindowSize)
		,	maxPairLzWindowSize(Default::MaxPairLzWindowSize)
		,	extraReduceHardReads(Default::ExtraReduceHardReads)
		,	extraReduceExpensiveLzMatches(Default::ExtraReduceExpensiveLzMatches)
	{}
};


struct ReadsContigBuilderParams
{
	struct Default
	{
		static const uint32 BeginCut = 2;
		static const uint32 EndCut = 2;
		static const uint32 MaxNewVariantsPerRead = 1;
		static const uint32 MaxRecordShiftDiff = 0;			// INFO: 0 - auto, half of the length (tests: 42)
		static const uint32 MaxHammingDistance = 8;			// TODO: a better metric?
		static const uint32 MinConsensusSize = 10;
	};

	uint32 beginCut;
	uint32 endCut;
	uint32 maxNewVariantsPerRead;
	uint32 maxRecordShiftDifference;
	uint32 maxHammingDistance;
	uint32 minConsensusSize;

	ReadsContigBuilderParams()
		:	beginCut(Default::BeginCut)
		,	endCut(Default::EndCut)
		,	maxNewVariantsPerRead(Default::MaxNewVariantsPerRead)
		,	maxRecordShiftDifference(Default::MaxRecordShiftDiff)
		,	maxHammingDistance(Default::MaxHammingDistance)
		,	minConsensusSize(Default::MinConsensusSize)
	{}
};


struct BinExtractorParams
{
	struct Default
	{
		static const uint32 MinBinSize = 256;				// 64 --> 512
	};

	uint32 minBinSize;

	BinExtractorParams()
		:	minBinSize(Default::MinBinSize)
	{}
};


struct CompressorParams
{
	struct Default
	{
		static const bool UseStoredToplogy = false;
		static const uint32 MaxMismatchesLowCost = 4;
	};


	MinimizerParameters minimizer;
	BinExtractorParams extractor;
	ReadsClassifierParams classifier;
	ReadsContigBuilderParams consensus;
	QualityCompressionParams quality;
	ArchiveType archType;

	bool useStoredTopology;
	uint32 maxMismatchesLowCost;

	CompressorParams()
		:	useStoredTopology(Default::UseStoredToplogy)
		,	maxMismatchesLowCost(Default::MaxMismatchesLowCost)
	{}
};


// INFO: now only for temporary file writing
#include <mutex>
#include <thread>

struct CompressorAuxParams
{
	bool dry_run;
	bool output_fastq;

	std::string uncompressed_filename;
	std::string uncompressed_filename_2;

	// Add a file to write reconstructed QV without the need of decompressing them
	FILE *f_uncompressed;
	FILE *f_uncompressed_2;

	// mutex for multithreaded writing in dry PE mode
	std::mutex* pe_mutex;

	CompressorAuxParams()
		:	dry_run(false)
		,	output_fastq(false)
		,	f_uncompressed(NULL)
		,	f_uncompressed_2(NULL)
		,	pe_mutex(NULL)
	{}
};


#endif // H_PACKPARAMS
