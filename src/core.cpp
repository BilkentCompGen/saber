#include "core.h"

seqan::String<seqan::DnaString> read_cores(seqan::CharString cores_filename) {   
    seqan::String<seqan::DnaString> core_list;
    
    std::ifstream infile(toCString(cores_filename));

    std::string core;
    while (std::getline(infile, core))
        appendValue(core_list, core);

    return core_list;    
}

int find_core_locations(pattern_locations &cores, seqan::DnaString &seq, aho_corasick_pattern &pattern) {
    int count = 0;
    seqan::Finder<seqan::DnaString> finder(seq);

    while (find(finder, pattern)) {
        int start = position(finder);
        int len = endPosition(finder) - position(finder) + 1;

        cores.emplace(std::make_pair(start, len));
        count++;
    }

    return count;
}

void print_core_locations(pattern_locations &cores) {
    std::cout << "{" << std::endl;
    for (std::pair<int,int> p : cores) {
       std::cout << " (" << p.first << ", " << p.second << ")" << std::endl;
    }
    std::cout << "}" << std::endl;    
}

bool is_core(std::pair<int,int> location, pattern_locations &cores) {
    return cores.find(location) != cores.end();
}
