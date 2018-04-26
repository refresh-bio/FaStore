/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef QVZ_H
#define QVZ_H

#include "../core/Globals.h"
#include "../core/BitMemory.h"

// QVZ includes
#include "Quality.h"		// quality compression params
#include "Stats.h"			// for QVZ stats
#include "pmf.h"
#include "well.h"
#include "distortion.h"
#include "quantizer.h"
#include "qv_file.h"
#include "codebook.h"

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
