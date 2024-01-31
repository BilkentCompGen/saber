#ifndef __dnablock_h
#define __dnablock_h

#include "types.h"
#include <stdexcept>

/**
 * Represents a block of a DNA sequence.
 */
class dna_block {
private:
    dna_sequence *seq; // Pointer to the DNA sequence this block belongs to.
    int start; // The start position of the block in the sequence (inclusive).
    int end; // The end position of the block in the sequence (exclusive).
    bool reversed; // Whether the block is reversed.

public:
    /**
     * Constructs a DNA block.
     * 
     * @param s Pointer to the DNA sequence this block belongs to.
     * @param start The start position of the block in the sequence.
     * @param end The end position of the block in the sequence.
     * @param r Whether the block is reversed (default is false).
     */
    dna_block(dna_sequence *s, int start, int end, bool r = false) 
        : seq(s), start(start), end(end), reversed(r) {
        if (s == nullptr) {
            throw std::invalid_argument("dna_sequence pointer cannot be nullptr");
        }
    }

    /**
     * Returns the size of the block.
     * 
     * @return The size of the block.
     */
    int size() {
        return end - start;
    }
    
    /**
     * Reverses the block.
     */
    void reverse() {
        reversed = !reversed;
    }

    /**
     * Checks whether the block is reversed.
     * 
     * @return True if the block is reversed, false otherwise.
     */
    bool is_reversed() {
        return reversed;
    }

    /**
     * Returns the range of the block in the sequence.
     * 
     * @return A pair of integers representing the start and end positions of the block.
     */
    std::pair<int, int> get_range() {
        return std::make_pair(start, end);
    }

    /**
     * Returns the character at the specified index in the block.
     * 
     * @param index The index of the character.
     * @return The character at the specified index.
     */
    char operator[](int index);

    /**
     * Returns the block as a DNA sequence.
     * 
     * @return The block as a DNA sequence.
     */
    dna_sequence get_block();

    /**
     * Returns the block as a C-style string.
     * 
     * @return The block as a C-style string.
     */
    char *c_str();
    
};

#endif //__dnablock_h