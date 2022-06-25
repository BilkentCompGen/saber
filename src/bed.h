#ifndef __bed__h
#define __bed__h

#include "core.h"
#include "edit_dist.h"

#include <algorithm>
#include <tuple>

// Approximate Block Edit Distance
int approximate_bed(seqan::DnaString &seq1, seqan::DnaString &seq2, aho_corasick_pattern &pattern);

// Direct Computation of Block Matches
void construct_match_table(std::vector<std::vector<int>> &d, dna_block &block, std::string &B);
void construct_match_table(std::vector<std::vector<int>> &d, std::string &_block, std::string &_B);

int find_starting_position(std::vector<std::vector<int>> &d, int _j, int best_end, dna_block &block, std::string &B);
int find_starting_position(std::vector<std::vector<int>> &d, int _j, int best_end, std::string &_block, std::string &_B);

// ignore matches criteria
bool ignore_match(int dist, int length1, int length2, double error_rate);
// Calculate W
void calculate_W(std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, seqan::DnaString &seq1, seqan::DnaString &seq2, int min_block, int max_block, double error_rate);
// Calculate M
void calculate_N(std::vector<int> &N, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, seqan::DnaString &seq1, seqan::DnaString &seq2, int max_iterations, int min_block, int max_block);

// Each block match includes the following information:
//      match[0] = start index of block from string A
//      match[1] = end index of block from string A (exclusive)
//      match[2] = start index of block from string B
//      match[3] = end index of block from string B (exclusive)
//      match[4] = 1 if reversed, else 0
//      match[5] = distance between blocks
//      
//      if indexes 2 and 3 are -1, then it is considered a block remove
void get_block_matches(std::vector<struct block_match> &block_matches, std::vector<int> &N, int i, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int min_block, int max_block);
void get_remaining_characters(std::string &new_seq1, std::string &new_seq2, seqan::DnaString &seq1, seqan::DnaString &seq2, int _i, std::vector<std::vector<int>> &block_matches);

int block_edit_score(std::vector<struct block_match> &block_matches, std::string &new_seq1, std::string &new_seq2, seqan::DnaString &seq1, seqan::DnaString &seq2, std::vector<int> &N, int _i, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int cutoff, int min_block, int max_block);
int block_edit_score(seqan::DnaString &seq1, seqan::DnaString &seq2, std::vector<int> &N, int i, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int cutoff, int min_block, int max_block);

int block_edit_score(seqan::DnaString &seq1, seqan::DnaString &seq2, std::vector<struct block_match> &block_matches);

bool block_overlaps_in_B(std::vector<int> &N, int _i, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int min_block, int max_block);

// overall block edit distance calculation.
int bed(std::vector<struct block_match> &block_matches, std::string &_seq1, std::string &_seq2, int max_iterations, int min_block, int max_block, double error_rate);

#endif // __bed__h