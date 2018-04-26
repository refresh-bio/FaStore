/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on QVZ software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef _CODEBOOK_H_
#define _CODEBOOK_H_
/**
 * Functions and definitions relating to reading codebooks from files, used
 * for both the encoder and decoder code
 */

#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <time.h>

#include "well.h"
#include "pmf.h"
#include "distortion.h"
#include "quantizer.h"
#include "qv_file.h"

#define MODE_RATIO		0	// Traditional implementation, output bitrate is scaled from input
#define MODE_FIXED		1	// Fixed rate per symbol
#define MODE_FIXED_MSE	2	// Fixed average MSE per column

#ifdef __MMAP__
#define qv2ch_mmap(a) a - 33
#else
#define qv2ch_mmap(a) a
#endif
#define qv2ch(a) a - 33

/**
 * Stores an array of conditional PMFs for the current column given the previous
 * column. PMF pointers are stored in a flat array so don't try to find the PMF you
 * want directly--use the accessor
 */
struct cond_pmf_list_t {
	uint32_t columns;
	const struct alphabet_t *alphabet;
	struct pmf_t **pmfs;
	struct pmf_list_t *marginal_pmfs;
    uint32_t pmfs_length;
};



// Memory management
//struct cond_pmf_list_t *alloc_conditional_pmf_list(const struct alphabet_t *alphabet, uint32_t columns);
struct cond_quantizer_list_t *alloc_conditional_quantizer_list(uint32_t columns);
void free_conditional_pmf_list(struct cond_pmf_list_t *);
void free_cond_quantizer_list(struct cond_quantizer_list_t *);

// Duplicate a cond_quant_list
struct cond_quantizer_list_t *duplicate_conditional_quantizer_list(struct cond_quantizer_list_t *in);
// Per-column initializer for conditional quantizer list
void cond_quantizer_init_column(struct cond_quantizer_list_t *list, uint32_t column, const struct alphabet_t *input_union);

// Accessors
struct pmf_t *get_cond_pmf(struct cond_pmf_list_t *list, uint32_t column, symbol_t prev);

uint32_t find_state_encoding(struct quantizer_t *codebook, symbol_t value);

void compute_qpmf_list(struct pmf_list_t *qpmf_list, struct cond_pmf_list_t *in_pmfs, uint32_t column, struct pmf_list_t *prev_qpmf_list, struct alphabet_t * q_alphabet_union, struct alphabet_t * prev_q_alphabet_union, struct cond_quantizer_list_t *q_list);
void compute_qpmf_quan_list(struct quantizer_t *q_lo, struct quantizer_t *q_hi, struct pmf_list_t *q_x_pmf, double ratio, struct alphabet_t *q_output_union);
void compute_xpmf_list(struct pmf_list_t *qpmf_list, struct cond_pmf_list_t *in_pmfs, uint32_t column, struct pmf_list_t *xpmf_list, struct alphabet_t * q_alphabet_union);
double optimize_for_distortion(struct pmf_t *pmf, struct distortion_t *dist, double target, struct quantizer_t **lo, struct quantizer_t **hi);

// Meat of the implementation
void calculate_statistics(struct quality_file_t *);
double optimize_for_entropy(struct pmf_t *pmf, struct distortion_t *dist, double target, struct quantizer_t **lo, struct quantizer_t **hi);
void generate_codebooks(struct quality_file_t *info);

// Master functions to handle codebooks in the output file
void write_codebooks(FILE *fp, struct quality_file_t *info);
void write_codebook(FILE *fp, struct cond_quantizer_list_t *quantizers);
void read_codebooks(FILE *fp, struct quality_file_t *info);
struct cond_quantizer_list_t *read_codebook(FILE *fp, struct quality_file_t *info);

#define MAX_CODEBOOK_LINE_LENGTH 3366
#define COPY_Q_TO_LINE(line, q, i, size) for (i = 0; i < size; ++i) { line[i] = q[i] + 33; }
#define COPY_Q_FROM_LINE(line, q, i, size) for (i = 0; i < size; ++i) { q[i] = line[i] - 33; }


void print_codebook(struct cond_quantizer_list_t *);

void preprocess_qvs(struct quality_file_t qv_info, struct qv_options_t *opts);
struct cluster_list_t *alloc_cluster_list(struct quality_file_t *info);

#endif
