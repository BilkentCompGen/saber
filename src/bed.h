#ifndef __bed__h
#define __bed__h

#include "ed.h"
#include "types.h"
#include "constants.h"

#include <algorithm>
#include <tuple>

/**
 * Constructs a match table for a block and a sequence. This is used to compute W
 * 
 * @param d The match table to be constructed.
 * @param block A block from source sequence.
 * @param B The target sequence.
 */
void construct_match_table(std::vector<std::vector<int>> &d, dna_block &block, dna_sequence &B);

/**
 * Backtracks through the match table to find the best starting for a match. This is used to compute W.
 * 
 * @param d The match table.
 * @param _j The current position in the sequence.
 * @param best_end The best ending position found so far.
 * @param block A block from source sequence.
 * @param B The target sequence.
 * @return The starting position of the match.
 */
int find_starting_position(std::vector<std::vector<int>> &d, int _j, int best_end, dna_block &block, dna_sequence &B);

/**
 * Determines whether a match should be ignored based on its distance and lengths.
 * 
 * @param dist The distance of the match.
 * @param l1 The length of the first block.
 * @param l2 The length of the second block.
 * @param error_rate The maximum error rate to be considered.
 * @return True if the match should be ignored, false otherwise.
 */
bool ignore_match(int dist, int l1, int l2, double error_rate);

/**
 * Calculates the W table for two sequences.
 * 
 * @param W The W table to be calculated.
 * @param seq1 The first sequence.
 * @param seq2 The second sequence.
 * @param settings The alignment settings.
 */
void calculate_W(std::vector<std::vector<w_entry>> &W, dna_sequence &seq1, dna_sequence &seq2, alignment_setting settings);

/**
 * Calculates the N table for two sequences.
 * 
 * @param N The N table to be calculated.
 * @param W The W table.
 * @param seq1 The first sequence.
 * @param seq2 The second sequence.
 * @param settings The alignment settings.
 */
void calculate_N(std::vector<int> &N, std::vector<std::vector<w_entry>> &W, dna_sequence &seq1, dna_sequence &seq2, alignment_setting settings);

/**
 * Gets the block matches for a sequence.
 * 
 * @param block_matches The block matches to be obtained.
 * @param N The N table.
 * @param i The current position in the sequence.
 * @param W The W table.
 * @param settings The alignment settings.
 */
void get_block_matches(std::vector<block_match> &block_matches, std::vector<int> &N, int i, std::vector<std::vector<w_entry>> &W, alignment_setting settings);

/**
 * Gets the remaining characters in two sequences after a match.
 * 
 * @param new_seq1 The remaining characters from first sequence to be calculated.
 * @param new_seq2 The remaining characters from second sequence to be calculated.
 * @param seq1 The original first sequence.
 * @param seq2 The original second sequence.
 * @param _i The current position in the sequences.
 * @param block_matches The block matches.
 */
void get_remaining_characters(dna_sequence &new_seq1, dna_sequence &new_seq2, dna_sequence &seq1, dna_sequence &seq2, int _i, std::vector<std::vector<int>> &block_matches);

/**
 * Calculates the block edit score for two sequences.
 * 
 * @param block_matches The block matches.
 * @param new_seq1 The first sequence after the match.
 * @param new_seq2 The second sequence after the match.
 * @param seq1 The original first sequence.
 * @param seq2 The original second sequence.
 * @param N The N table.
 * @param _i The current position in the sequences.
 * @param W The W table.
 * @param cutoff The cutoff score.
 * @param settings The alignment settings.
 * @return The block edit score.
 */
int block_edit_score(std::vector<block_match> &block_matches, dna_sequence &new_seq1, dna_sequence &new_seq2, dna_sequence &seq1, dna_sequence &seq2, std::vector<int> &N, int _i, std::vector<std::vector<w_entry>> &W, int cutoff, alignment_setting settings);

/**
 * Calculates the block edit score for two sequences.
 * 
 * @param seq1 The first sequence.
 * @param seq2 The second sequence.
 * @param N The N table.
 * @param i The current position in the sequences.
 * @param W The W table.
 * @param cutoff The cutoff score.
 * @param settings The alignment settings.
 * @return The block edit score.
 */
int block_edit_score(dna_sequence &seq1, dna_sequence &seq2, std::vector<int> &N, int i, std::vector<std::vector<w_entry>> &W, int cutoff, alignment_setting settings);

/**
 * Calculates the block edit score for two sequences.
 * 
 * @param seq1 The first sequence.
 * @param seq2 The second sequence.
 * @param block_matches The block matches.
 * @return The block edit score.
 */
int block_edit_score(dna_sequence &seq1, dna_sequence &seq2, std::vector<block_match> &block_matches);

/**
 * Determines whether a block overlaps in sequence B.
 * 
 * @param N The N table.
 * @param _i The current position in the sequence.
 * @param W The W table.
 * @param settings The alignment settings.
 * @return True if the block overlaps, false otherwise.
 */
bool block_overlaps_in_B(std::vector<int> &N, int _i, std::vector<std::vector<w_entry>> &W, alignment_setting settings);

/**
 * Calculates the overall block edit distance for two sequences.
 * 
 * @param block_matches The block matches to be computed.
 * @param _seq1 The first sequence.
 * @param _seq2 The second sequence.
 * @param max_iterations The maximum number of iterations.
 * @param settings The alignment settings.
 * @return The block edit distance.
 */
int bed(std::vector<block_match> &block_matches, dna_sequence &_seq1, dna_sequence &_seq2, int max_iterations, alignment_setting settings);

#endif // __bed__h