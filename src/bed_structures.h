#ifndef __bed_structures_h
#define __bed_structures_h

#include <iostream>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>

struct block_match {
    std::pair<int,int> loc1;
    std::pair<int,int> loc2;
    bool remove;
    bool reversed;
    int dist;
};

// Blocks Information
// TODO: MAKE THESE THE PARAMETERS
//#define MIN_BLOCK_LEN 50 // should be at least 2
//#define MAX_BLOCK_LEN 75
#define BLOCK_RANGE (MAX_BLOCK_LEN - MIN_BLOCK_LEN + 1)

#define CONSIDER_MOVE_CONSTANT 0.0 // positive. recommended: (0.1)
//#define ERROR_RATE 0.3 // in (0, 1). recommended: 0.2

// Unknown distance between between blocks
#define INF_DIST INT_MAX
#define UNKNOWN_TUPLE std::make_tuple(-1,-1, false)

// Character Operation Costs
#define C_CHAR_DELETION 1
#define C_CHAR_INSERTION 1
#define C_CHAR_SNP 1

// Block Operation Costs
#define C_BLOCK_DELETION 1
#define C_BLOCK_MOVE 1
#define C_BLOCK_REVERSE 1

// typedef for aho_corasick_pattern
typedef seqan::Pattern<seqan::String<seqan::DnaString>, seqan::AhoCorasick> aho_corasick_pattern; 
typedef std::set<std::pair<int, int>> pattern_locations;

// typedef for a dna block
typedef seqan::Infix<seqan::DnaString>::Type dna_block;

#endif // __bed_structures_h