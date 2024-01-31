#ifndef __ed_h
#define __ed_h

#include "types.h"
#include "dnablock.h"

#include "edlib/edlib/include/edlib.h"

// Levenshtein dist from edlib.
/**
 * Computes the Levenshtein distance between two strings using the Edlib library.
 * 
 * @param _a First string.
 * @param m Length of the first string.
 * @param _b Second string.
 * @param n Length of the second string.
 * @param k Maximum allowed distance. If the actual distance is greater than k, the function returns -1.
 * @return The Levenshtein distance between _a and _b, or -1 if the distance is greater than k.
 */
int lev_dist_edlib(const char *_a, int m, const char *_b, int n, int k);

/**
 * Computes the Levenshtein distance between two dna_block objects using the Edlib library.
 * 
 * @param _a First dna_block object.
 * @param _b Second dna_block object.
 * @param k Maximum allowed distance. If the actual distance is greater than k, the function returns -1.
 * @return The Levenshtein distance between the sequences represented by _a and _b, or -1 if the distance is greater than k.
 */
int lev_dist_edlib(dna_block &_a, dna_block &_b, int k);

/**
 * Computes the Levenshtein distance between two dna_sequence objects using the Edlib library.
 * 
 * @param _a First dna_sequence object.
 * @param _b Second dna_sequence object.
 * @param k Maximum allowed distance. If the actual distance is greater than k, the function returns -1.
 * @return The Levenshtein distance between the sequences represented by _a and _b, or -1 if the distance is greater than k.
 */
int lev_dist_edlib(dna_sequence &_a, dna_sequence &_b, int k);

#endif // __ed_h