/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_PACKPARAMS
#define H_PACKPARAMS

#include "../core/Globals.h"

#include "ContigBuilder.h"

#include "../core/FastqCategorizer.h"
#include "../core/ReadsClassifier.h"
#include "../qvz/Quality.h"
#include "../fastore_bin/Params.h"
#include "../fastore_rebin/BinFileExtractor.h"


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
