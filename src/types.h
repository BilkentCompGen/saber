#ifndef __types_h
#define __types_h

#include <utility>
#include <string>

class block_match {
public:
    std::pair<int,int> loc1;
    std::pair<int,int> loc2;
    bool remove;
    bool reversed;
    int dist;
};

typedef std::string dna_sequence;

#endif // __types_h