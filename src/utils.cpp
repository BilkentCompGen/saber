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
    {'N', 'N'},
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
    if (c == 'N') return base_lookup[dis(gen)];
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
    std::ifstream file(filename, std::ios::binary);
    if (!file)
        return false;

    // Get file size to pre-allocate memory
    file.seekg(0, std::ios::end);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    sequence.clear();
    if (size > 0) {
        sequence.reserve(size);
    }

    std::string line;
    enum FileType { AUTO, FASTA, FASTQ, RAW } type = AUTO;

    while (std::getline(file, line)) {
        // Remove trailing carriage return, newline, or spaces
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' ')) {
            line.pop_back();
        }

        // Find first non-whitespace character
        size_t first_non_space = line.find_first_not_of(" \t\r\n");
        if (first_non_space == std::string::npos) {
            continue; // Skip empty lines
        }
        if (first_non_space > 0) {
            line = line.substr(first_non_space);
        }

        if (type == AUTO) {
            if (line[0] == '>') {
                type = FASTA;
                continue;
            } else if (line[0] == '@') {
                type = FASTQ;
                continue;
            } else {
                type = RAW;
            }
        }

        if (type == FASTA) {
            if (line[0] == '>') {
                // Hit next sequence header, stop reading
                break;
            }
        }

        // Convert sequence to uppercase and clean internal whitespace
        std::transform(line.begin(), line.end(), line.begin(), ::toupper);
        line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char c) {
            return std::isspace(c);
        }), line.end());

        // Convert non-ATGCN characters to 'N'
        for (char &c : line) {
            if (c != 'A' && c != 'T' && c != 'G' && c != 'C' && c != 'N') {
                c = 'N';
            }
        }

        sequence.append(line);

        if (type == FASTQ) {
            // For FASTQ, we only read the single sequence line and then stop
            break;
        }
    }

    file.close();
    return !sequence.empty();
}
