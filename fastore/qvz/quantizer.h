/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on QVZ software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef _QUANTIZER_H_
#define _QUANTIZER_H_

#include <stdint.h>

#include "pmf.h"
#include "distortion.h"
#include "util.h"
#include "well.h"

#define QUANTIZER_MAX_ITER		100

/**
 * Structure holding information about a quantizer, which just maps input symbols
 * to output symbols for a specific alphabet
 */
struct quantizer_t {
	const struct alphabet_t *alphabet;
	struct alphabet_t *output_alphabet;
	symbol_t *q;
    double ratio;
	double mse;
};

/**
 * Stores an array of quantizer pointers for the column for all possible left context
 * values. Unused ones are left as null pointers. This is also stored as a flat array
 * so the accessor must be used to look up the correct quantizer
 * The dreaded triple pointer is used to store an array of (different length) arrays
 * of pointers to quantizers
 */
struct cond_quantizer_list_t {
    uint32_t columns;
    uint32_t lines;
    struct alphabet_t **input_alphabets;
    struct quantizer_t ***q;
    double **ratio;				// Raw ratio
    uint8_t **qratio;			// Quantized ratio
    const struct qv_options_t *options;
};

// Memory management
struct quantizer_t *alloc_quantizer(const struct alphabet_t *);
void free_quantizer(struct quantizer_t *);

// Generates a quantizer via optimization
struct quantizer_t *generate_quantizer(struct pmf_t *pmf, struct distortion_t *dist, uint32_t states);

// Calculate the output pmf when the quantizer is applied to the input pmf
struct pmf_t *apply_quantizer(struct quantizer_t *q, struct pmf_t *pmf, struct pmf_t *output);

// Find the output alphabet of a quantizer
void find_output_alphabet(struct quantizer_t *);

// Accessors
struct quantizer_t *get_cond_quantizer_indexed(struct cond_quantizer_list_t *list, uint32_t column, uint32_t index);
struct quantizer_t *get_cond_quantizer(struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev);
void store_cond_quantizers(struct quantizer_t *lo, struct quantizer_t *hi, double ratio, struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev);
void store_cond_quantizers_indexed(struct quantizer_t *lo, struct quantizer_t *hi, double ratio, struct cond_quantizer_list_t *list, uint32_t column, uint32_t index);
struct quantizer_t *choose_quantizer(struct cond_quantizer_list_t *list, struct well_state_t *well, uint32_t column, symbol_t prev, uint32_t *q_idx);

// Display/debugging
void print_quantizer(struct quantizer_t *);

#endif
