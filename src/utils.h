#ifndef __utils_h
#define __utils_h

#include <random>
#include "types.h"

char complement(char c);

class coin_tosser {
private:
    std::mt19937 gen;
    std::uniform_int_distribution<> dis;

public:
    // Constructor that seeds the random number generator with the current time
    coin_tosser() : gen(std::random_device{}()), dis(0, 1) {}

    // Perform a coin toss
    bool toss() {
        return dis(gen) == 0;
    }
};


class random_base_generator {
private:
    std::mt19937 gen;
    std::uniform_int_distribution<> dis;

public:
    // Constructor that seeds the random number generator with the current time
    random_base_generator() : gen(std::random_device{}()), dis(0, 4) {}

    // Generate base
    char generate();
};

class base_mutator {
private:
    std::mt19937 gen;
    std::uniform_int_distribution<> dis;

public:
    // Constructor that seeds the random number generator with the current time
    base_mutator() : gen(std::random_device{}()), dis(0, 3) {}

    // Generate base
    char mutate(char c);
};  

bool read_sequence(const std::string& filename, dna_sequence& sequence);

#endif //__utils_h