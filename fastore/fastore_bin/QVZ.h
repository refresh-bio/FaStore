/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef QVZ_H
#define QVZ_H

#include "../fastore_bin/Globals.h"
#include "../fastore_bin/Quality.h"		// quality compression params
#include "../fastore_bin/Stats.h"			// for QVZ stats

#include "BitMemory.h"

// QVZ includes
#include "../fastore_pack/pmf.h"
#include "../fastore_pack/well.h"
#include "../fastore_pack/distortion.h"
#include "../fastore_pack/quantizer.h"
#include "../fastore_pack/qv_file.h"
#include "../fastore_pack/codebook.h"

struct cond_quantizer_list_t;
struct qv_options_t;


// QVZ-codebook calculation
//
struct QvzCodebook
{
	struct cond_quantizer_list_t *qlist;

	QvzCodebook();
	~QvzCodebook();

	void ComputeFromStats(cond_pmf_list_t* trainingStats_, const struct qv_options_t *qvzOpts);
	void WriteCodebook(BitMemoryWriter& fp, uint32 max_columns);
	void ReadCodebook(BitMemoryReader& fp, struct alphabet_t *in_alphabet, uint32_t columns);
};



// essential data for quality compression used when compressing
// and decompressing the q-scores
// (apart from quality compression parameters)
//
struct QualityCompressionData
{
	QvzCodebook codebook;
    well_state_t well;
	uint32 max_read_length;
};



#endif // QVZ_H
