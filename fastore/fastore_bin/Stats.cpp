/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "Stats.h"

// QVZ includes
#include "../fastore_pack/pmf.h"
#include "../fastore_pack/well.h"
#include "../fastore_pack/distortion.h"
#include "../fastore_pack/quantizer.h"
#include "../fastore_pack/qv_file.h"
#include "../fastore_pack/codebook.h"


FastqRawBlockStats::FastqRawBlockStats()
{
	qua.columns =  FastqRecord::MaxSeqLen;
	qua.training_stats = alloc_conditional_pmf_list(ALPHABET_SIZE, qua.columns);

	Clear();

}

FastqRawBlockStats::~FastqRawBlockStats()
{
	// TODO: cleaning of training_stats
}

void FastqRawBlockStats::Clear()
{
	std::fill(dna.symFreq.begin(), dna.symFreq.end(), 0);
	std::fill(qua.symFreq.begin(), qua.symFreq.end(), 0);

	// We need only to clear total and counts because we have still not
	// not compute the propoer pmfs or marginal_pmfs
	uint32_t alphabet_card = qua.training_stats->alphabet->size;
	uint32_t pmfs_length = qua.training_stats->pmfs_length;

	for(uint32_t i = 0; i < pmfs_length; ++i)
	{
		qua.training_stats->pmfs[i]->total = 0;
		ASSERT(qua.training_stats->pmfs[i]->pmf_ready == 0);

		for(uint32_t j = 0; j < alphabet_card; ++j )
		{
			qua.training_stats->pmfs[i]->counts[j] = 0;
		}
	}


	// clear headers
	//
	head.fields.clear();
}

// updates stats per-record while processing records
//
void FastqRawBlockStats::Update(const FastqRecord &rec_)
{
	FastqRecordBinStats::Update(rec_);

	// HINT: for clarity this can be also moved into
	// DnaStats / QualityStats :: Update() fcn
	//
	for (uint32 i = 0; i < rec_.seqLen; ++i)
		dna.symFreq[rec_.seq[i]]++;


	if (rec_.qua != NULL)
	{
		pmf_increment(get_cond_pmf(qua.training_stats, 0, 0), qv2ch(rec_.qua[0]));
		qua.symFreq[rec_.qua[0]]++;

		for (uint32 i = 1; i < rec_.seqLen; ++i)
		{
			qua.symFreq[rec_.qua[i]]++;
			pmf_increment(get_cond_pmf(qua.training_stats, i, qv2ch(rec_.qua[i-1])), qv2ch(rec_.qua[i]) );
		}
	}



	// are we using headers?
	//
	if (rec_.head != NULL)
	{
		ASSERT(rec_.headLen > 0);

		// check tokens -- TODO: optimize
		//
		const std::string separators = FastqRawBlockStats::HeaderStats::Separators();

		uint32 fieldNo = 0;
		uint32 fieldStartPos = 0;

		for (uint32 i = 0; i <= rec_.headLen; ++i)
		{
			if (!std::count(separators.begin(), separators.end(), rec_.head[i]) && (i != rec_.headLen))
				continue;

			// check whether the field has been already set
			//
			const char* field = rec_.head + fieldStartPos;
			const uint32 fieldLen = i - fieldStartPos;
			if (head.fields.size() < fieldNo + 1)
			{
				// setup new field
				//
				head.fields.push_back(HeaderStats::Field());
				HeaderStats::Field& f = head.fields.back();
				f.isConst = true;

				uint64 v;
				if (is_num(field, fieldLen, v))
				{
					f.isNumeric = true;
					f.minValue = f.maxValue = v;
				}
				else
				{
					f.isNumeric = false;
					f.minValue = f.maxValue = fieldLen;
					f.possibleValues.insert(std::string(field, fieldLen));
				}


				// add new separator
				//
				if (i != rec_.headLen)
					f.separator = rec_.head[i];
			}
			else
			{
				// update the existing field
				//
				uint64 v;
				bool isNumeric = is_num(field, fieldLen, v);

				HeaderStats::Field& f = head.fields[fieldNo];
				ASSERT(f.isNumeric == isNumeric);

				if (isNumeric)
				{
					f.minValue = std::min(f.minValue, v);
					f.maxValue = std::max(f.maxValue, v);
					f.isConst &= (f.minValue == f.maxValue);
				}
				else
				{
					f.possibleValues.insert(std::string(field, fieldLen));
					f.isConst &= f.possibleValues.size() == 1;
				}


				// verify the separator
				//
				if (i != rec_.headLen)
					ASSERT(f.separator == rec_.head[i]);
			}

			fieldStartPos = i + 1;
			fieldNo++;
		}
	}
}


// updates stats after processing bins
//
void FastqRawBlockStats::Update(const FastqRawBlockStats &stats_)
{
	FastqRecordBinStats::Update(stats_);

	// HINT: for clarity this can be also moved into
	// DnaStats / QualityStats :: Update() fcn
	//
	for (uint32 i = 0; i < stats_.dna.symFreq.size(); ++i)
		dna.symFreq[i] += stats_.dna.symFreq[i];


	// update quality
	//
	for (uint32 i = 0; i < stats_.qua.symFreq.size(); ++i)
		qua.symFreq[i] += stats_.qua.symFreq[i];


	uint32_t alphabet_card = stats_.qua.training_stats->alphabet->size;
	uint32_t pmfs_length = stats_.qua.training_stats->pmfs_length;

	for(uint32_t i = 0; i < pmfs_length; ++i)
	{
		qua.training_stats->pmfs[i]->total += stats_.qua.training_stats->pmfs[i]->total;
		for(uint32_t j = 0; j < alphabet_card; ++j )
		{
			qua.training_stats->pmfs[i]->counts[j] += stats_.qua.training_stats->pmfs[i]->counts[j];
		}
	}



	// update headers
	//
	if (stats_.head.fields.size() > 0)
	{
		ASSERT(stats_.head.fields.size() == head.fields.size() || head.fields.size() == 0);
		if (head.fields.size() == 0)
		{
			head.fields = stats_.head.fields;
		}
		else
		{
			for (uint32 i = 0; i < stats_.head.fields.size(); ++i)
			{
				HeaderStats::Field& f1 = head.fields[i];
				const HeaderStats::Field& f2 = stats_.head.fields[i];

				ASSERT(f1.isNumeric == f2.isNumeric);
				ASSERT(f1.separator == f2.separator);
				if (f1.isNumeric)
				{
					f1.minValue = std::min(f1.minValue, f2.minValue);
					f1.maxValue = std::max(f1.maxValue, f2.maxValue);
					f1.isConst &= f1.minValue == f1.maxValue;
				}
				else
				{
					f1.possibleValues.insert(f2.possibleValues.begin(), f2.possibleValues.end());
					f1.isConst &= f1.possibleValues.size() == 1;
				}
			}
		}
	}
}


struct cond_pmf_list_t * FastqRawBlockStats::alloc_conditional_pmf_list(uint32_t alphabet_size, uint32_t columns)
{
	uint32_t count = 1 + alphabet_size*(columns-1);
	uint32_t i;
	struct cond_pmf_list_t *list = (struct cond_pmf_list_t *) calloc(1, sizeof(struct cond_pmf_list_t));

	// We need one array of PMF pointers that will index into the buffer allocated above, for the columns
	list->columns = columns;
	list->alphabet = alloc_alphabet(alphabet_size);
	list->pmfs = (struct pmf_t **) calloc(count, sizeof(struct pmf_t *));
	list->pmfs_length = count;
	// All PMFs are stored in a flat array, the accessor function will resolve a PMF's address
	for (i = 0; i < count; ++i) {
		list->pmfs[i] = alloc_pmf(list->alphabet);
	}

	return list;
}


void FastqRawBlockStats::Compute_marginal_pmf()
{
	cond_pmf_list_t *pmf_list = qua.training_stats;

	pmf_list->marginal_pmfs = alloc_pmf_list(qua.columns, pmf_list->alphabet);
	combine_pmfs(get_cond_pmf(pmf_list, 0, 0), pmf_list->marginal_pmfs->pmfs[0], 1.0, 0.0, pmf_list->marginal_pmfs->pmfs[0]);
	for (uint32_t column = 1; column < qua.columns; ++column) {
		for (uint32_t j = 0; j < pmf_list->alphabet->size; ++j) {
			combine_pmfs(pmf_list->marginal_pmfs->pmfs[column], get_cond_pmf(pmf_list, column, j), 1.0, get_probability(pmf_list->marginal_pmfs->pmfs[column-1], j), pmf_list->marginal_pmfs->pmfs[column]);
		}
	}
}
