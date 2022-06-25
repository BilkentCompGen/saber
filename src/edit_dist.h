#ifndef __edit_dist_h
#define __edit_dist_h

#include "bed_structures.h"
#include "edlib/edlib/include/edlib.h"

template <class T>
T min(const T& left, const T& right) {
	return left > right ? right : left;
}

template <class T>
T max(const T& left, const T& right) {
	return left > right ? left : right;
}

// levenshtein edit distance implementation
int lev_dist(seqan::DnaString &_a, seqan::DnaString &_b);
int lev_dist(dna_block &_a, dna_block &_b);
int lev_dist(dna_block &_a, std::string &_b);
int lev_dist(std::string &_a, std::string &_b);
int lev_dist(dna_block &_a, seqan::DnaString &_b);

// Levenshtein edit distance implementation using Myers' Bit Vector Algorithm, 
// only when both strings are less than 64 characters long.
int lev_dist_myers_64(std::string &_a, std::string &_b);

// Levenshtein dist from edlib.
int lev_dist_edlib(std::string &_a, std::string &_b, int k);
int lev_dist_edlib(seqan::DnaString &_a, seqan::DnaString &_b, int k);


// reverse complement operation implementation
std::string reverse_complement(seqan::DnaString &seq);
std::string reverse_complement(dna_block &seq);
std::string reverse_complement(std::string &seq);

// Ttansform to std string from seqan strings
std::string to_string(seqan::DnaString &seq);
std::string to_string(dna_block &seq);

// edit distance including reverse complement operation
int dist_rev(bool &_reversed, seqan::DnaString &_a, seqan::DnaString &_b);

int dist_rev(bool &_reversed, dna_block &_a, dna_block &_b);
int dist_rev(bool &_reversed, dna_block &_a, dna_block &_b, std::string &_b_rev);

int dist_rev(bool &_reversed, std::string &_a, std::string &_b);


#endif // __edit_dist_h