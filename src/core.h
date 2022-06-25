// Ahmet Cemal Alıcıoğlu

#ifndef __core_h
#define __core_h

// STD includes
#include <string>
#include <fstream>
#include <unordered_set>
#include <sys/time.h>
#include <ctime>

// SeqAn includes
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>

// Local includes
#include "bed_structures.h"


// struct for hashing of core locations
/*
struct core_location_hash {
    inline std::size_t operator()(const std::pair<int,int> &v) const {
        return 1000 * (v.second - MIN_CORE_LEN) + v.first;
    }
}*/

seqan::String<seqan::DnaString> read_cores(seqan::CharString cores_filename);
int find_core_locations(pattern_locations &cores, seqan::DnaString &seq, aho_corasick_pattern &pattern);
void print_core_locations(pattern_locations &cores);
bool is_core(std::pair<int,int> location, pattern_locations &cores);

#endif // __core_h