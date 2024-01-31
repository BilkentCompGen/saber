#ifndef __types_h
#define __types_h

#include "constants.h"

#include <utility>
#include <string>

/**
 * Represents a match between two blocks of sequences.
 */
class block_match {
public:
    std::pair<int,int> loc1;  // Location of the first block (start, end)
    std::pair<int,int> loc2;  // Location of the second block (start, end)
    bool remove;  // Whether this is a block remove or not (loc2 and reversed is ignored if true)
    bool reversed;  // Whether the match is a reverse
    int dist;  // The distance between the two blocks


    block_match(int loc1_start, int loc1_end) : loc1(std::make_pair(loc1_start, loc1_end)), loc2(std::make_pair(-1,-1)), remove(true), reversed(false), dist(0) {}
    block_match(int loc1_start, int loc1_end, int loc2_start, int loc2_end, bool _reversed, int _dist) : loc1(std::make_pair(loc1_start, loc1_end)), loc2(std::make_pair(loc2_start, loc2_end)), remove(false), reversed(_reversed), dist(_dist) {}
};


/**
 * Represents an entry in a dynamic programming table for sequence alignment.
 */
class w_entry {
    int score;  // The score of this entry
    int start;  // The start position of the block in the target sequence
    int end;  // The end position of the block in the target sequence
    bool reversed;  // Whether the block is reversed

public:
    int get_score() { return score; }
    int get_start() { return start; }
    int get_end() { return end; }
    bool is_reversed() { return reversed; }
    std::pair<int,int> get_block() { return std::make_pair(start, end); }

    w_entry() : score(INF_DIST), start(0), end(0), reversed(false) {}    
    w_entry(int _score, int _start, int _end, bool _reversed) : score(_score), start(_start), end(_end), reversed(_reversed) {}

    void set(int _score, int _start, int _end, bool _reversed) {
        score = _score;
        start = _start;
        end = _end;
        reversed = _reversed;
    }

    bool valid() { return score != INF_DIST; }
};

/**
 * Represents the settings for a sequence alignment.
 */
class alignment_setting {
public:
    int min_block;  // The minimum size of a block for the alignment
    int max_block;  // The maximum size of a block for the alignment
    int max_iterations;  // The maximum number of iterations for the alignment algorithm
    double error_rate;  // The maximum character error rate considered for the block matches

    alignment_setting(int _min_block = DEFAULT_L_MIN, int _max_block = DEFAULT_L_MAX, int _max_iterations = DEFAULT_MAX_ITERATIONS, double _error_rate = DEFAULT_ERROR_RATE) : min_block(_min_block), max_block(_max_block), max_iterations(_max_iterations), error_rate(_error_rate) {}
};

typedef std::string dna_sequence;

#endif // __types_h