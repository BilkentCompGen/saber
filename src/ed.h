#ifndef __ed_h
#define __ed_h

#include "types.h"
#include "dnablock.h"

#include "edlib/edlib/include/edlib.h"

// Levenshtein dist from edlib.
int lev_dist_edlib(const char *_a, int m, const char *_b, int n, int k);
int lev_dist_edlib(dna_block &_a, dna_block &_b, int k);
int lev_dist_edlib(dna_sequence &_a, dna_sequence &_b, int k);

#endif // __ed_h