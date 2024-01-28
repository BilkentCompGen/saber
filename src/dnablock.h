#ifndef __dnablock_h
#define __dnablock_h

#include "types.h"
#include <stdexcept>


class dna_block {
private:
    dna_sequence *seq;
    int start; // inclusive
    int end; // not inclusive
    bool reversed;

public:
    dna_block(dna_sequence *s, int start, int end, bool r = false) 
        : seq(s), start(start), end(end), reversed(r) {
        if (s == nullptr) {
            throw std::invalid_argument("dna_sequence pointer cannot be nullptr");
        }
    }

    int size() {
        return end - start;
    }
    
    void reverse() {
        reversed = !reversed;
    }

    bool is_reversed() {
        return reversed;
    }

    std::pair<int, int> get_range() {
        return std::pair<int,int>(start, end);
    }

    char operator[](int index);
    dna_sequence get_block();
    char *c_str();
    
};

#endif //__dnablock_h