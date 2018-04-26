/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on QVZ software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef _WELL_H_
#define _WELL_H_
/**
 * Declarations necessary for using the WELL-1024a PRNG in other code
 */

#include <stdint.h>

struct well_state_t {
	uint32_t state[32];
	uint32_t n;
	uint32_t bit_output;
	uint32_t bits_left;
};


uint32_t well_1024a(struct well_state_t *state);
uint32_t well_1024a_bits(struct well_state_t *state, uint8_t bits);

#endif
