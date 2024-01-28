#include "dnablock.h"
#include "utils.h"


char dna_block::operator[](int index) {
    if (index < 0 || index >= size()) 
        throw std::out_of_range("Index out of range");
                
    if (reversed) {
        return complement((*seq)[end - index - 1]);
    } else {
        return (*seq)[start + index];
    }
}


dna_sequence dna_block::get_block() {
    int n = size();
    dna_sequence result(n, 0);

    for (int i = 0; i < n; i++) 
        result[i] = (*this)[i];
        
    return result;
}

char *dna_block::c_str() {
    int n = size();
    char *str = new char[n + 1];
        
    for (int i = 0; i < n; i++) 
        str[i] = (*this)[i];
        
    str[n] = '\0';  // null-terminate the string
        
    return str;
}