#ifndef __utils_h
#define __utils_h

#include <random>
#include "types.h"

/**
 * Returns the complement of a DNA base.
 * 
 * @param c The DNA base.
 * @return The complement of the DNA base.
*/
char complement(char c);

/**
 * Coin tosser class
 */
class coin_tosser {
private:
    std::mt19937 gen;
    std::uniform_int_distribution<> dis;

public:
    // Constructor that seeds the random number generator with the current time
    coin_tosser() : gen(std::random_device{}()), dis(0, 1) {}

    /**
     * Toss a coin
     * @return true or false randomly
    */
    bool toss() {
        return dis(gen) == 0;
    }
};

/**
 * Random DNA base generator 
*/
class random_base_generator {
private:
    std::mt19937 gen;
    std::uniform_int_distribution<> dis;

public:
    // Constructor that seeds the random number generator with the current time
    random_base_generator() : gen(std::random_device{}()), dis(0, 3) {}

    /**
     * Generate base
     * @return random base (A, T, C, G)
    */
    char generate();
};

/**
 * SNP simulator
*/
class base_mutator {
private:
    std::mt19937 gen;
    std::uniform_int_distribution<> dis;

public:
    // Constructor that seeds the random number generator with the current time
    base_mutator() : gen(std::random_device{}()), dis(1, 3) {}

    /**
     * Mutates a DNA base
     * @param c the base
     * @return mutated base
    */
    char mutate(char c);
};  

/**
 * Reads a DNA sequence from a file.
 * 
 * @param filename The name of the file to read from.
 * @param sequence The DNA sequence to read into.
 * @return True if the sequence was successfully read, false otherwise.
 */
bool read_sequence(const std::string& filename, dna_sequence& sequence);

/**
 * Returns the minimum of three values.
 * 
 * @param a The first value.
 * @param b The second value.
 * @param c The third value.
 * @return The minimum of the three values.
 */
template <typename T>
T min_3(const T& a, const T& b, const T& c) {
    return std::min(a, std::min(b, c));
}

#endif //__utils_h