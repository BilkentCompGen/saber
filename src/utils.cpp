#include "utils.h"

#include <stdexcept>
#include <fstream>
#include <map>
#include <iostream>
#include <algorithm>
#include <cctype>

static const std::map<char, char> complement_lookup = {
    {'A', 'T'},
    {'T', 'A'},
    {'C', 'G'},
    {'G', 'C'},
};

static const char base_lookup[4] = {'A', 'T', 'C', 'G'};

static const std::map<char, int> base_lookup_rev = {
    {'A', 0},
    {'T', 1},
    {'C', 2},
    {'G', 3},
};

char random_base_generator::generate() {
    return base_lookup[dis(gen)];
}

char base_mutator::mutate(char c) {
    if (base_lookup_rev.find(c) == base_lookup_rev.end())
        throw std::invalid_argument("Invalid nucleotide");

    return base_lookup[(base_lookup_rev.at(c) + dis(gen)) % 4];
}

char complement(char c) {
    if (complement_lookup.find(c) == complement_lookup.end())
        throw std::invalid_argument("Invalid nucleotide");
    
    return complement_lookup.at(c);
}

bool read_sequence(const std::string& filename, dna_sequence& sequence) {
    std::ifstream file(filename);
    if (!file)
        return false;

    std::string line;
    while (std::getline(file, line)) {
        // Skip the header line in FASTA/FASTQ
        if (line.empty() || line[0] == '>' || line[0] == '@') continue;

        // Read the sequence line
        sequence = line;
        std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        sequence.erase(std::remove_if(sequence.begin(), sequence.end(), [](char c) {
            return c == '\n' || c == ' ';
        }), sequence.end());

        // Check for invalid characters
        if (std::any_of(sequence.begin(), sequence.end(), [](char c) {
            return c != 'A' && c != 'T' && c != 'G' && c != 'C';
        })) {
            file.close();
            return false;
        }

        break;
    }

    file.close();
    return true;
}
