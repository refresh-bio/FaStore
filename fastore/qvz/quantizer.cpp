/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on QVZ software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "quantizer.h"
#include "util.h"

/**
 *
 */
int compare_ptr_double (const void * a, const void * b)
{
    double ** v1 = (double**)a;
    double ** v2 = (double**)b;
    if (**v2 > **v1)
        return 1;
    else if (**v1 > **v2)
        return -1;
    else
        return 0;
    //return (int)( **v2 - **v1 );
}
/**
 *
 */
int compare_int (const void * a, const void * b)
{
    symbol_t* v1 = (symbol_t*)a;
    symbol_t* v2 = (symbol_t*)b;
    if (*v2 > *v1)
        return -1;
    else if (*v1 > *v2)
        return 1;
    else
        return 0;
    //return (int)( **v2 - **v1 );
}
/**
 * Allocate enough room based on the size of the alphabet supplied
 */
struct quantizer_t *alloc_quantizer(const struct alphabet_t *alphabet) {
    struct quantizer_t *rtn = (struct quantizer_t *) calloc(1, sizeof(struct quantizer_t));
	rtn->alphabet = alphabet;
	rtn->q = (symbol_t *) calloc(alphabet->size, sizeof(symbol_t));
	return rtn;
}

/**
 * Free the quantizer itself but not the input alphabet
 * But do free the output alphabet
 */
void free_quantizer(struct quantizer_t *q) {
	if (q->output_alphabet)
		free_alphabet(q->output_alphabet);

	free(q->q);
	free(q);
}

/**
 * Produce a quantizer with the given number of states for the given pmf, and
 * optionally computes the expected distortion produced by this quantizer.
 * The bounds array here contains the left endpoint (inclusive) of each region
 */
struct quantizer_t *generate_quantizer(struct pmf_t *pmf, struct distortion_t *dist, uint32_t states) {
	struct quantizer_t *q = alloc_quantizer(pmf->alphabet);
	uint32_t changed = 1;
	uint32_t iter = 0;
	int32_t i, j, r, size;
	uint32_t min_r;
	double mse, min_mse, next_mse;
	symbol_t *bounds = (symbol_t *) _alloca((states+1)*sizeof(symbol_t));
	symbol_t *reconstruction = (symbol_t *) _alloca(states*sizeof(symbol_t));

    // OLD VERSION
	// Initial bounds and reconstruction points
    /*
	bounds[0] = 0;
	bounds[states] = pmf->alphabet->size;
	for (j = 1; j < states; ++j) {
		bounds[j] = (j * pmf->alphabet->size) / states;
	}
	for (j = 0; j < states; ++j) {
		reconstruction[j] = (bounds[j] + bounds[j+1] - 1) / 2;
	}
    */
    size = pmf->alphabet->size;
    int mass_count = 0;
    int first_mass = -1, last_mass = 0;
    int *mass_array = (int*)calloc(size, sizeof(int));
    // Vesrion 4.2
    for (i=0; i<size; i++) {
        if (pmf->pmf[i] > 0)
        {
            if(first_mass < 0){
                first_mass = i;
                
            }
            mass_array[mass_count] = i;
            mass_count++;
            last_mass = i;
        }
    }
    if (mass_count == 0 || states == 1) {
        assert(states==1);
        bounds[0] = 0;
        bounds[states] = ALPHABET_SIZE;
    }
    
    else{
        bounds[states] = last_mass+1;
        for (j=0; j<states; j++) {
            bounds[j] = mass_array[(j * mass_count) / states];
        }
        
    }
    for (j = 0; j < states; ++j) {
        reconstruction[j] = (bounds[j] + bounds[j+1] - 1) / 2;
    }
    free(mass_array);
    /*
    // Version 4.1
    // Initial bounds and reconstruction points
    int mass_count = 0;
    int first_mass = -1, last_mass = 0;
    for (i=0; i<size; i++) {
        if (pmf->pmf[i] > 0)
        {
            if(first_mass < 0){
                first_mass = i;
                
            }
            mass_count++;
            last_mass = i;
        }
    }
    if (mass_count == 0 && states == 1) {
        bounds[0] = 0;
        bounds[states] = ALPHABET_SIZE;
    }
    else if (mass_count == 0 && states > 1){
        // This should not happend
        printf("");
    }
    else{
        bounds[0] = first_mass;
        bounds[states] = last_mass + 1;
        size = (last_mass + 1) - first_mass;
        for (j = 1; j < states; ++j) {
            bounds[j] = ((j * size) / states)+first_mass;
        }
    }
    for (j = 0; j < states; ++j) {
        reconstruction[j] = (bounds[j] + bounds[j+1] - 1) / 2;
    }
    
    size = pmf->alphabet->size;
    
    */
    
    /* Version 3
    int mass_count = 0;
    bounds[0] = 0;
    bounds[states] = pmf->alphabet->size;
    for (j = 1; j < states; ++j) {
        bounds[j] = (j * pmf->alphabet->size) / states;
    }
    while (1) {
        
        mass_count = 0;
        for (i=bounds[j-1]; i<bounds[j]; i++) {
            if (pmf->pmf[i] > 0)
                mass_count++;
        }
        while ((mass_count += pmf->pmf[i++]) == 0) {
            bounds[j]++;
        }
        
    }
    for (j = 0; j < states; ++j) {
        reconstruction[j] = (bounds[j] + bounds[j+1] - 1) / 2;
    }
    
    size = pmf->alphabet->size;
    */
    /*
    double target_mass = 1.0/states;
    int set_count = 0;
    double cum_mass = 0.0;
    int max_count, mass_count=0, mass_count_ctr = 0;
    
    for (i=0; i<size; i++) {
        if (pmf->pmf[i] > 0)
            mass_count++;
    }
    if (mass_count == 0) {
        // I think this should not happend.
        // We need to go through the code to try to figure it out
        printf("");
    }
    
    if (mass_count == 0 && states > 1) {
        // I think this should not happend.
        // We need to go through the code to try to figure it out
        printf("");
    }
    
    max_count = mass_count - states + 1;
    
    //assert(max_count > 0);
    bounds[0] = 0;
    bounds[states] = pmf->alphabet->size;
    j=0;
    for (i=0; i<size && mass_count_ctr <= mass_count; i++) {
        cum_mass += pmf->pmf[i];
        if (pmf->pmf[i] > 0) {
            set_count++;
            mass_count_ctr++;
            if (set_count == max_count) {
                bounds[++j] = i+1;
                cum_mass = 0.0;
                set_count = 0;
                max_count = mass_count - mass_count_ctr - states + j + 1;
                continue;
            }
        }
        if (cum_mass >= target_mass) {
            bounds[++j] = i+1;
            max_count = mass_count - mass_count_ctr - states + j + 1;
            cum_mass = 0.0;
            set_count = 0;
        }
    }
    // First, adjust the reconstruction points for fixed bounds
    for (j = 0; j < states; ++j) {
        // Initial guess for min values
        min_mse = DBL_MAX;
        min_r = bounds[j];
        assert(bounds[j] < ALPHABET_SIZE);
        assert(bounds[j+1] <= ALPHABET_SIZE);
        // For each possible reconstruction point
        for (r = bounds[j]; r < bounds[j+1]; ++r) {
            // Find its distortion when used for the whole region
            mse = 0.0;
            for (i = bounds[j]; i < bounds[j+1]; ++i) {
                mse += get_probability(pmf, i) * get_distortion(dist, i, r);
            }
            
            // Compare to minimums, save if better
            if (mse < min_mse) {
                min_r = r;
                min_mse = mse;
            }
        }
        
        // Check if we've changed our reconstruction and save it
        if (min_r != reconstruction[j]) {
            changed = 1;
            reconstruction[j] = min_r;
        }
        assert(reconstruction[j] < ALPHABET_SIZE);
    }
    */
    // VERSION 2
    /*
    // This is a k-means algorithm, thus initialization is very important.
    // We are setting the initial points to those with the largest mass.
    // First we sort the and keep track of the indexes
    
    // Create a pointer array to the pmfs
    double** pmf_ptr = (double**)_alloca(size*sizeof(double*));
    //double *pmf_ptr[ALPHABET_SIZE];
    for (i = 0; i < size; i++) {
        pmf_ptr[i] = &(pmf->pmf[i]);
    }
    // Sort the pointer array. The index of the sorted values are (pmf_ptr[i] - pmf->pmf)
    qsort(pmf_ptr, size, sizeof(*pmf_ptr), compare_ptr_double);
    // Assign a reconstruction point to the most probable symbols
    for (j = 0; j < states; ++j) {
        reconstruction[j] = pmf_ptr[j] - pmf->pmf;
        //reconstruction[j] = j;
    }
    // Sort again the reconstruction vector.
    qsort(reconstruction, states, sizeof(*reconstruction), compare_int);
    //free(pmf_ptr);
    // Then, adjust the bounds for fixed reconstruction points by iterating
    // over the positions (apart from the endpoints which always have a fixed
    // assignment) and deciding which of the two nearest points they
    // contribute the least expected distortion to
    bounds[0] = 0;
    bounds[states] = size;
    r = 0;
    for (j = 1; j < size && r < states-1; ++j) {
        // Get distortion for the current and next reconstruction points
        // I don't think the PMF actually affects this since it is the same
        // coefficient for both and we are comparing them
        mse = get_distortion(dist, j, reconstruction[r]);
        next_mse = get_distortion(dist, j, reconstruction[r+1]);
        
        // if the next one is lower, save the current symbol as the left bound
        // for that region
        if (next_mse < mse) {
            r += 1;
            bounds[r] = j;
            assert(bounds[r] < ALPHABET_SIZE);
        }
    }
    for (j = 0; j <= states; ++j) {
        assert(bounds[j] <= ALPHABET_SIZE);
    }
    */
    
	// Lloyd-Max quantizer design alternating between adjustment of bounds
	// and of reconstruction point locations until there is no change
	while (changed && iter < QUANTIZER_MAX_ITER) {
		changed = 0;
		iter += 1;

		// First, adjust the reconstruction points for fixed bounds
		for (j = 0; j < states; ++j) {
			// Initial guess for min values
			min_mse = DBL_MAX;
			min_r = bounds[j];
            assert(bounds[j] < ALPHABET_SIZE);
            assert(bounds[j+1] <= ALPHABET_SIZE);
			// For each possible reconstruction point
			for (r = bounds[j]; r < bounds[j+1]; ++r) {
				// Find its distortion when used for the whole region
				mse = 0.0;
				for (i = bounds[j]; i < bounds[j+1]; ++i) {
					mse += get_probability(pmf, i) * get_distortion(dist, i, r);
				}

				// Compare to minimums, save if better
				if (mse < min_mse) {
					min_r = r;
					min_mse = mse;
				}
			}

			// Check if we've changed our reconstruction and save it
			if (min_r != reconstruction[j]) {
				changed = 1;
				reconstruction[j] = min_r;
			}
            assert(reconstruction[j] < ALPHABET_SIZE);
		}

		// Then, adjust the bounds for fixed reconstruction points by iterating
		// over the positions (apart from the endpoints which always have a fixed
		// assignment) and deciding which of the two nearest points they
		// contribute the least expected distortion to
		r = 0;
		for (j = 1; j < size-1 && r < states-1; ++j) {
			// Get distortion for the current and next reconstruction points
			// I don't think the PMF actually affects this since it is the same
			// coefficient for both and we are comparing them
			mse = get_distortion(dist, j, reconstruction[r]);
			next_mse = get_distortion(dist, j, reconstruction[r+1]);

			// if the next one is lower, save the current symbol as the left bound
			// for that region
			if (next_mse < mse) {
				r += 1;
				bounds[r] = j;
			}
		}
	}
    
    // We need to initialize the quantizer to NOT_SYMBOLS for the decoder
    for (i = 0; i<q->alphabet->size; i++) {
        q->q[i] = ALPHABET_NOT_SYMBOL;
    }

	// Now, iterate over the regions and set up the quantizer mapping from input
	// to reconstruction point
	for (j = 0; j < states; ++j) {
		for (i = bounds[j]; i < bounds[j+1]; ++i) {
			q->q[i] = reconstruction[j];
		}
	}
    
    // Add this to make the initialization of the Lloyd-max invisible to the decoder
    for (i = 0; i<bounds[0]; i++) {
        q->q[i] = q->q[bounds[0]];
    }
    for (i = bounds[states]; i<q->alphabet->size; i++) {
        q->q[i] = q->q[bounds[states-1]];
    }

	// Save the output alphabet in the quantizer
	q->output_alphabet = alloc_alphabet(states);
	memcpy(q->output_alphabet->symbols, reconstruction, sizeof(symbol_t) * states);
	alphabet_compute_index(q->output_alphabet);

	// Calculate the distortion and store it in the quantizer
	q->mse = 0.0;
	for (j = 0; j < states; ++j) {
		for (i = bounds[j]; i < bounds[j+1]; ++i) {
			q->mse += get_distortion(dist, i, reconstruction[j]) * get_probability(pmf, i);
		}
	}
    
	return q;
}

/**
 * Calculate the PMF of the output when the given quantizer is used with symbols generated
 * from the given input distribution. The input and output pmf structures cannot be the
 * same. If output is null, a new PMF will be allocated and a pointer returned
 */
struct pmf_t *apply_quantizer(struct quantizer_t *q, struct pmf_t *pmf, struct pmf_t *output) {
	uint32_t i;

	if (!pmf->pmf_ready)
		recalculate_pmf(pmf);
	
	if (output) {
		// Clear existing pmf from output
		memset(output->pmf, 0, output->alphabet->size * sizeof(double));
	}
	else {
		// Allocate a new PMF for output
		output = alloc_pmf(pmf->alphabet);
	}

	// Sum together input probabilities that map to the same output
	for (i = 0; i < pmf->alphabet->size; ++i) {
		output->pmf[q->q[i]] += get_probability(pmf, i);
	}
	output->pmf_ready = 1;

	return output;
}

/**
 * Generates the output alphabet from the quantization table, in case this isn't
 * already available
 */
void find_output_alphabet(struct quantizer_t *q) {
	symbol_t p;
	uint32_t x;
	uint32_t size = 0;
	symbol_t *uniques = (symbol_t *) _alloca(q->alphabet->size * sizeof(symbol_t));

	// First symbol in quantizer output is always unique
	p = q->q[0];
    if (p != ALPHABET_NOT_SYMBOL) {
        uniques[0] = p;
        size = 1;
    }
	

	// Search the rest of the quantizer
    // We needed to change this since now the end of the quantizer may be ALPHABET_NOT_SYMBOLs (unused values)
	for (x = 1; x < q->alphabet->size; ++x) {
		if (q->q[x] != p) {
			p = q->q[x];
            if (p == ALPHABET_NOT_SYMBOL) {
                break;
            }
			uniques[size] = p;
			size += 1;
		}
	}

	// Make it into a proper alphabet
	q->output_alphabet = alloc_alphabet(size);
	memcpy(q->output_alphabet->symbols, uniques, size*sizeof(symbol_t));
	alphabet_compute_index(q->output_alphabet);
}

/**
 * Get a quantizer by its indexed location within the quantizer list for a column
 */
struct quantizer_t *get_cond_quantizer_indexed(struct cond_quantizer_list_t *list, uint32_t column, uint32_t index) {
    return list->q[column][index];
}

/**
 * Get a quantizer by its left context symbol
 */
struct quantizer_t *get_cond_quantizer(struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev) {
    uint32_t idx = get_symbol_index(list->input_alphabets[column], prev);
    if (idx != ALPHABET_SYMBOL_NOT_FOUND)
        return get_cond_quantizer_indexed(list, column, idx);
    return NULL;
}

/**
 * Stores the given quantizers at the appropriate index corresponding to the left context symbol given
 * for the specific column
 */
void store_cond_quantizers(struct quantizer_t *lo, struct quantizer_t *hi, double ratio, struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev) {
    uint32_t idx = get_symbol_index(list->input_alphabets[column], prev);
    store_cond_quantizers_indexed(lo, hi, ratio, list, column, idx);
}

/**
 * Stores the given quantizers directly at the given index based on the previous context symbol. Faster when
 * we know what the previous index was in addition to what the previous symbol was
 */
void store_cond_quantizers_indexed(struct quantizer_t *lo, struct quantizer_t *hi, double ratio, struct cond_quantizer_list_t *list, uint32_t column, uint32_t idx) {
    list->q[column][2*idx] = lo;
    list->q[column][2*idx + 1] = hi;
    list->ratio[column][idx] = ratio;
    list->qratio[column][idx] = (uint8_t) (ratio * 128.);
}

/**
 * Selects a quantizer for the given column from the quantizer list with the appropriate ratio
 */
struct quantizer_t *choose_quantizer(struct cond_quantizer_list_t *list, struct well_state_t *well, uint32_t column, symbol_t prev, uint32_t *q_idx) {
    uint32_t idx = get_symbol_index(list->input_alphabets[column], prev);
    assert(idx != ALPHABET_SYMBOL_NOT_FOUND);
    if (well_1024a_bits(well, 7) >= list->qratio[column][idx]) {
        *q_idx = 2*idx+1;
        return list->q[column][2*idx+1];
    }
    *q_idx = 2*idx;
    return list->q[column][2*idx];
}


/**
 * Print a quantizer to stdout
 */
void print_quantizer(struct quantizer_t *q) {
	uint32_t i;
	char *tmp = (char *) _alloca(q->alphabet->size+1);

	tmp[q->alphabet->size] = 0;
	for (i = 0; i < q->alphabet->size; ++i) {
		tmp[i] = (char) (q->q[i] + 33);
	}
	printf("Quantizer: %s\n", tmp);

	tmp[q->output_alphabet->size] = 0;
	for (i = 0; i < q->output_alphabet->size; ++i) {
		tmp[i] = (char) (q->output_alphabet->symbols[i] + 33);
	}
	printf("Unique alphabet: %s\n", tmp);
}
