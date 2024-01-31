#include "dnablock.h"
#include "utils.h"


char dna_block::operator[](int index) {
    if (index < 0 || index >= size()) 
        throw std::out_of_range("Index out of range");
                
    if (reversed)
        return complement((*seq)[end - index - 1]);
     
    return (*seq)[start + index];
    
}

dna_sequence dna_block::get_block() {
    int n = size();
    dna_sequence result;
    result.reserve(n);

    if (reversed) {
        for (int i = 0; i < n; i++) {
            result[i] = complement((*seq)[end - i - 1]);
        }
    } else {
        for (int i = 0; i < n; i++) {
            result[i] = (*seq)[start + i];
        }
    }

    return result;
}

char *dna_block::c_str() {
    int n = size();
    char *str = new char[n + 1];

    if (reversed) {
        for (int i = 0; i < n; i++) {
            str[i] = complement((*seq)[end - i - 1]);
        }
    } else {
        for (int i = 0; i < n; i++) {
            str[i] = (*seq)[start + i];
        }
    }

    str[n] = '\0';  // null-terminate the string

    return str;
}