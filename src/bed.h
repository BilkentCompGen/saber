#ifndef __bed__h
#define __bed__h

#include "ed.h"
#include "types.h"

#include <algorithm>
#include <tuple>


// Direct Computation of Block Matches
void construct_match_table(std::vector<std::vector<int>> &d, dna_block &block, dna_sequence &B);

int find_starting_position(std::vector<std::vector<int>> &d, int _j, int best_end, dna_block &block, dna_sequence &B);

// ignore matches criteria
bool ignore_match(int dist, int l1, int l2, double error_rate);

// Calculate W
void calculate_W(std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, dna_sequence &seq1, dna_sequence &seq2, int min_block, int max_block, double error_rate);
// Calculate M
void calculate_N(std::vector<int> &N, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, dna_sequence &seq1, dna_sequence &seq2, int max_iterations, int min_block, int max_block);

void get_block_matches(std::vector<struct block_match> &block_matches, std::vector<int> &N, int i, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int min_block, int max_block);
void get_remaining_characters(std::string &new_seq1, std::string &new_seq2, dna_sequence &seq1, dna_sequence &seq2, int _i, std::vector<std::vector<int>> &block_matches);

int block_edit_score(std::vector<struct block_match> &block_matches, std::string &new_seq1, std::string &new_seq2, dna_sequence &seq1, dna_sequence &seq2, std::vector<int> &N, int _i, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int cutoff, int min_block, int max_block);
int block_edit_score(dna_sequence &seq1, dna_sequence &seq2, std::vector<int> &N, int i, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int cutoff, int min_block, int max_block);

int block_edit_score(dna_sequence &seq1, dna_sequence &seq2, std::vector<struct block_match> &block_matches);

bool block_overlaps_in_B(std::vector<int> &N, int _i, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int min_block, int max_block);

// overall block edit distance calculation.
int bed(std::vector<struct block_match> &block_matches, std::string &_seq1, std::string &_seq2, int max_iterations, int min_block, int max_block, double error_rate);

#endif // __bed__h