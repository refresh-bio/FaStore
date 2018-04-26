/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on QVZ software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "QVZ.h"

#include <memory>

QvzCodebook::QvzCodebook()
	:	qlist(NULL)
{}

QvzCodebook::~QvzCodebook()
{
	if (qlist != NULL)
	{
		free_cond_quantizer_list(qlist);
	}
}


void QvzCodebook::ComputeFromStats(cond_pmf_list_t* trainingStats_, const struct qv_options_t *qvzOpts)
{
	// Stuff for state allocation and mixing
	double ratio;

	// Miscellaneous variables
	uint32_t column, j;
	double total_mse;
	double target_dist;

	uint32_t columns = trainingStats_->columns;

	// Output list of conditional quantizers
	struct cond_quantizer_list_t *q_list;

	// Constant alphabet of all possible input symbols
	const struct alphabet_t *A = trainingStats_->alphabet;

	// Temporary/extra pointers
	struct quantizer_t *q_lo;
	struct quantizer_t *q_hi;

	// List of conditionally quantized PMFs after quantizer has been added out
	struct pmf_list_t *xpmf_list;

	// List of conditionally quantized PMFs after the next quantizer was applied
	struct pmf_list_t *qpmf_list;
	struct pmf_list_t *prev_qpmf_list;

	// Alphabet of all possible quantizer outputs from the previous column
	struct alphabet_t *q_output_union;
	struct alphabet_t *q_prev_output_union;

	struct cond_pmf_list_t *in_pmfs;
	const struct qv_options_t *opts = qvzOpts;

	// Compute the distortion matrix
	struct distortion_t *dist;
	if (opts->distortion == DISTORTION_CUSTOM) {
		dist = gen_custom_distortion(ALPHABET_SIZE, opts->dist_file);
	}
	else {
		dist = generate_distortion_matrix(ALPHABET_SIZE, opts->distortion);
	}
	target_dist = opts->D;

	q_list = alloc_conditional_quantizer_list(columns);
	qlist = q_list;
	in_pmfs = trainingStats_;

	// For the column 0 the quantizers aren't conditional, so find them directly
	q_output_union = alloc_alphabet(1);
	cond_quantizer_init_column(q_list, 0, q_output_union);
	q_list->options = opts;

	// Initialize the new pmfs (dummy)
	qpmf_list = alloc_pmf_list(A->size, q_output_union);

	// Handle column zero specially
	// @todo handle fixed mse target
	//if (opts->mode == MODE_RATIO)
	//	ratio = optimize_for_entropy(get_cond_pmf(in_pmfs, 0, 0), dist, get_entropy(get_cond_pmf(in_pmfs, 0, 0))*opts->ratio, &q_lo, &q_hi);
	//else{
	//h1 = compute_DR_slope(get_cond_pmf(in_pmfs, 0, 0), dist);
	ratio = optimize_for_distortion(get_cond_pmf(in_pmfs, 0, 0), dist, opts->D, &q_lo, &q_hi);
	//}
	q_lo->ratio = ratio;
	q_hi->ratio = 1-ratio;
	total_mse = ratio*q_lo->mse + (1-ratio)*q_hi->mse;
	store_cond_quantizers(q_lo, q_hi, ratio, q_list, 0, 0);

	// free the used pmfs and alphabet
	// (do not free q_prev_output_union and prev_qpmf_output as it's the first assignment).
	q_prev_output_union = q_output_union;
	prev_qpmf_list = qpmf_list;

	// Start computing the quantizers of the rest of the columns
	for (column = 1; column < columns; column++) {
		// Compute the next output alphabet union over all quantizers for this column
		q_output_union = duplicate_alphabet(get_cond_quantizer_indexed(q_list, column-1, 0)->output_alphabet);
		for (j = 1; j < 2*q_prev_output_union->size; ++j) {
			alphabet_union(q_output_union, get_cond_quantizer_indexed(q_list, column-1, j)->output_alphabet, q_output_union);
		}
		cond_quantizer_init_column(q_list, column, q_output_union);

		// Initialize the new pmfs
		qpmf_list = alloc_pmf_list(A->size, q_output_union);
		xpmf_list = alloc_pmf_list(q_output_union->size, A);

		// Compute P(Q_i|X_i)
		if (column == 1)
			compute_qpmf_quan_list(q_lo, q_hi, qpmf_list, ratio, q_output_union);
		else
			compute_qpmf_list(qpmf_list, in_pmfs, column, prev_qpmf_list, q_output_union, q_prev_output_union, q_list);

		// Compute P(X_{i+1}|Q_i)
		compute_xpmf_list(qpmf_list, in_pmfs, column, xpmf_list, q_output_union);

		// for each previous value Q_i compute the quantizers
		for (j = 0; j < q_output_union->size; ++j) {
			// Find and save quantizers
			// @todo handle fixed mse target
			//if (opts->mode == MODE_RATIO)
			//	ratio = optimize_for_entropy(xpmf_list->pmfs[j], dist, get_entropy(xpmf_list->pmfs[j])*opts->ratio, &q_lo, &q_hi);
			//else{
			// For version 2
			//hi = compute_DR_slope(xpmf_list->pmfs[j], dist);
			//if (hi>0) {
			//    hi = h1;
			//}
			//target_dist = (hi/h1)*opts->D;
			///////
			target_dist = opts->D;
			ratio = optimize_for_distortion(xpmf_list->pmfs[j], dist,target_dist , &q_lo, &q_hi);
			//}
			q_lo->ratio = ratio;
			q_hi->ratio = 1-ratio;
			store_cond_quantizers_indexed(q_lo, q_hi, ratio, q_list, column, j);

			// This actually needs to be scaled by the probability of this quantizer pair being used to be accurate, uniform assumption is an approximation
			total_mse = (ratio*q_lo->mse + (1-ratio)*q_hi->mse);
		}

		// deallocated the memory of the used pmfs and alphabet
		free(q_prev_output_union);
		q_prev_output_union = q_output_union;
		free_pmf_list(prev_qpmf_list);
		prev_qpmf_list = qpmf_list;
		free_pmf_list(xpmf_list);
	}

	// Final cleanup, things we saved at the end of the final iteration that aren't needed
	free_pmf_list(qpmf_list);
	free(q_output_union);
	//}
}


void QvzCodebook::WriteCodebook(BitMemoryWriter& fp, uint32 max_columns)
{
	struct cond_quantizer_list_t* quantizers = qlist;
	uint32_t i, j, k;

	ASSERT(max_columns <= quantizers->columns);
	//uint32_t columns = quantizers->columns;
	uint32_t columns = max_columns;

	struct quantizer_t *q_temp = get_cond_quantizer_indexed(quantizers, 0, 0);
	uint32_t size = q_temp->alphabet->size;
	uint32_t buflen = columns > size ? columns : size;
	//char *eol = "\n";
	char *linebuf = (char *) _alloca(sizeof(char)*buflen);

	// First chunk, ratio for zero context quantizer
	fp.PutByte(quantizers->qratio[0][0] + 33);

	// Second chunk is low quantizer
	fp.Put2Bytes(size);				// store the size so we won't need to use EOLs

	COPY_Q_TO_LINE(linebuf, q_temp->q, i, size);
	fp.PutBytes((byte*)linebuf, sizeof(char) * size);

	// Third line is high quantizer
	q_temp = get_cond_quantizer_indexed(quantizers, 0, 1);
	COPY_Q_TO_LINE(linebuf, q_temp->q, i, size);
	fp.PutBytes((byte*)linebuf, sizeof(char) * size);


	// Now for the rest of the columns, use the same format
	for (i = 1; i < columns; ++i) {
		// First a line containing ratios for each previous context
		for (j = 0; j < quantizers->input_alphabets[i]->size; ++j) {
			linebuf[j] = quantizers->qratio[i][j] + 33;
		}
		// store the alphabet size so we won't need to use EOLs
		fp.Put2Bytes(quantizers->input_alphabets[i]->size);
                //ASSERT(quantizers->input_alphabets[i]->size == size);

		fp.PutBytes((byte*)linebuf, sizeof(char) * quantizers->input_alphabets[i]->size);

		// Next, the low quantizers in index order
		for (j = 0; j < quantizers->input_alphabets[i]->size; ++j) {
			q_temp = get_cond_quantizer_indexed(quantizers, i, 2*j);
			COPY_Q_TO_LINE(linebuf, q_temp->q, k, size);
			fp.PutBytes((byte*)linebuf, sizeof(char) * size);
		}

		// Finally, the high quantizers in index order
		for (j = 0; j < quantizers->input_alphabets[i]->size; ++j) {
			q_temp = get_cond_quantizer_indexed(quantizers, i, 2*j+1);
			COPY_Q_TO_LINE(linebuf, q_temp->q, k, size);
			fp.PutBytes((byte*)linebuf, sizeof(char) * size);
		}
	}
}


void QvzCodebook::ReadCodebook(BitMemoryReader& fp, struct alphabet_t *in_alphabet, uint32_t columns) {
   uint32_t column, size, out_alphabet_size, part_size;
   uint32_t i, j;
   struct quantizer_t *q_lo, *q_hi;
   struct alphabet_t *uniques;
   char line[MAX_CODEBOOK_LINE_LENGTH];
   uint8_t qratio;
   struct alphabet_t *A = in_alphabet;

   uniques = alloc_alphabet(1);
   qlist = alloc_conditional_quantizer_list(columns);
   cond_quantizer_init_column(qlist, 0, uniques);
   free_alphabet(uniques);

   // Next chunk is qratio for zero quantizer offset by 33
   qratio = fp.GetByte() - 33;

   out_alphabet_size = fp.Get2Bytes();
   ASSERT(out_alphabet_size == A->size);

   // Allocate some quantizers and copy the tables from lines 3 and 4
   q_lo = alloc_quantizer(A);
   q_hi = alloc_quantizer(A);
   fp.GetBytes((byte*)line, out_alphabet_size);
   COPY_Q_FROM_LINE(line, q_lo->q, j, A->size);
   fp.GetBytes((byte*)line, out_alphabet_size);
   COPY_Q_FROM_LINE(line, q_hi->q, j, A->size);

   // Fill in missing uniques information and store
   find_output_alphabet(q_lo);
   find_output_alphabet(q_hi);
   uniques = alloc_alphabet(0);
   alphabet_union(q_lo->output_alphabet, q_hi->output_alphabet, uniques);
   store_cond_quantizers_indexed(q_lo, q_hi, 0.0, qlist, 0, 0);
   qlist->qratio[0][0] = qratio;

   // Now handle the remaining columns uniformly
   for (column = 1; column < columns; ++column) {
	   // Initialize the column information so we can write to it directly
	   cond_quantizer_init_column(qlist, column, uniques);
	   size = uniques->size;
	   free_alphabet(uniques);
	   uniques = alloc_alphabet(0);

	   part_size = fp.Get2Bytes();

	   // First line is the ratios
	   fp.GetBytes((byte*)line, part_size);

	   for (i = 0; i < size; ++i) {
		   qlist->qratio[column][i] = line[i] - 33;
	   }

	   // Next line is a number of low quantizers
	   for (i = 0; i < size; ++i) {
		   q_lo = alloc_quantizer(A);
		   fp.GetBytes((byte*)line, A->size*sizeof(symbol_t));
		   COPY_Q_FROM_LINE(line, q_lo->q, j, A->size);

		   find_output_alphabet(q_lo);
		   qlist->q[column][2*i] = q_lo;
		   alphabet_union(uniques, q_lo->output_alphabet, uniques);
	   }

	   // Next line is a number of high quantizers
	   for (i = 0; i < size; ++i) {
		   q_hi = alloc_quantizer(A);
		   fp.GetBytes((byte*)line, A->size*sizeof(symbol_t));
		   COPY_Q_FROM_LINE(line, q_hi->q, j, A->size);

		   find_output_alphabet(q_hi);
		   qlist->q[column][2*i+1] = q_hi;
		   alphabet_union(uniques, q_hi->output_alphabet, uniques);
	   }
   }

   // We don't use the uniques from the last column
   free_alphabet(uniques);
}
