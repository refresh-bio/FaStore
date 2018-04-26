/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_STATS
#define H_STATS

#include "Globals.h"

#include <vector>
#include <set>
#include <algorithm>

#include "FastqRecord.h"
#include "Buffer.h"
#include "Utils.h"


/**
 * FASTQ records statistics which will be gathered during parsing
 *
 */
struct FastqRawBlockStats : public FastqRecordBinStats
{
	struct DnaStats
	{
		std::array<uint64, 128> symFreq;			// not used, just a placeholder
	};

	struct QualityStats
	{
		std::array<uint64, 128> symFreq;			// not used, just a placeholder

		struct cond_pmf_list_t *training_stats;		// used by QVZ to calculate codebooks
        uint32_t columns;
        struct qv_options_t *opts;
	};

	struct HeaderStats
	{
		static const uint32 MaxPossibleValues = 1 << 16;

		static const std::string Separators()
		{
			return " ./:#+";
		}

		struct Field
		{
			bool isConst;
			bool isNumeric;
			char separator;
			uint64 minValue;						// or min/max length
			uint64 maxValue;
			std::set<std::string> possibleValues;

			Field()
				:	isConst(false)
				,	isNumeric(false)
				,	separator(0)
				,	minValue((uint64)-1)
				,	maxValue(0)
			{}
		};

		std::vector<Field> fields;

		uint32 pairedEndFieldIdx;


		HeaderStats()
			:	pairedEndFieldIdx(0)
		{}
	};


	static const uint32 MaxSeqLen = FastqRecord::MaxSeqLen;

	DnaStats dna;
	QualityStats qua;
	HeaderStats head;


	FastqRawBlockStats();
	~FastqRawBlockStats();

	void Clear();
    
	struct cond_pmf_list_t * alloc_conditional_pmf_list(uint32_t alphabet_size, uint32_t columns);

	// updates stats per-record while processing records
	//
	void Update(const FastqRecord &rec_);

	// updates stats after processing bins
	//
	void Update(const FastqRawBlockStats &stats_);

	void Compute_marginal_pmf();
    
};

#endif
