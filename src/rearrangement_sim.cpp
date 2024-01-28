#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <random>
#include <string>
#include <fstream>

#include "ed.h"
#include "bed.h"
#include "constants.h"
#include "utils.h"

int mutate_seq(dna_sequence &mutated_seq, dna_sequence seq, double reverse_rate, double max_char_error_rate, int l_min, int l_max, int max_removes, int max_moves) {
    random_base_generator rbg;
    base_mutator bm;
    
    // Sequence Rearrangement Simulator
    int l_range = l_max - l_min + 1;
    int m = seq.length();

    std::srand(std::time(nullptr));

    int max_blocks = max_removes + max_moves;
    if (max_blocks == 0) {
        mutated_seq = seq;
        return 0;
    }
    
    std::vector<std::pair<int, int>> blocks;

    int range_per_block = m / max_blocks;

    std::string out_char = "";

    for (int i = 0; i < max_blocks; i++) {
        int len =  l_min + std::rand() % l_range;
        int max_end = range_per_block - len;
        int start = std::rand() %  max_end;
        blocks.emplace_back(range_per_block * i + start, range_per_block * i + start + len - 1);
    }

    for (int i = 0; i < m; i++) {
        bool out = true;
        for (int j = 0; j < max_blocks; j++) {
            if (i >= blocks[j].first && i <= blocks[j].second) {
                out = false;
                break;
            }
        }
        if (out)
            out_char.push_back(seq[i]);
    }

    int reverses = 0;
    std::vector<std::pair<int,bool>> v(max_blocks); // (move location, reversed?)
    for (int i = 0; i < max_moves; i++) {
        v[i].first = std::rand() % out_char.size();
        int num = std::rand() % 1000;
        if (num < reverse_rate * 1000) {
            v[i].second = true;
            reverses++;
        } else {
            v[i].second = false;
        }
    }

    for (int i = max_moves; i < max_blocks; i++) {
        v[i].first = -1;
    }

    std::random_shuffle(v.begin(), v.end());

    std::string new_seq = "";
    std::vector<struct block_match> block_matches;
    for (int i = 0; i < (int) out_char.size(); i++) {
        for (int j = 0; j < max_blocks; j++) {
            if (v[j].first == i) {
                std::cout << "\nm: ";
                if (!v[j].second) {
                    for (int k = blocks[j].first; k <= blocks[j].second; k++) {
                        new_seq.push_back(seq[k]);
                        std::cout << seq[k];
                    }

                } else {
                    for (int k = blocks[j].second; k >= blocks[j].first; k--) {
                        new_seq.push_back(complement(seq[k]));
                    }

                }
                std::cout << std::endl;
            }
        }
        new_seq.push_back(out_char[i]);
        std::cout << out_char[i];
    }
    std::cout << "\n" << std::endl;
    std::cout << "seq:\t\t" << seq << std::endl;
    std::cout << "new_seq:\t" << new_seq << std::endl;

    mutated_seq = "";
    int n = new_seq.size();

    int char_edits = 0;
    // Then char edits
    for (int i = 0; i < n; i++) {
        
        int num = std::rand() % 1000;
        if (num < (max_char_error_rate / 3) * 1000) {
            // snp
            mutated_seq.push_back(bm.mutate(new_seq[i]));

            char_edits++;
        } else if ((max_char_error_rate / 3) * 1000 <= num && num < (max_char_error_rate / 3) * 2000) {
            // in
            mutated_seq.push_back(rbg.generate());
            mutated_seq.push_back(new_seq[i]);

            char_edits++;
        } else if ((max_char_error_rate / 3) * 2000 <= num && num < (max_char_error_rate) * 1000) {
            // del
            char_edits++;
        } else {
            // no error
            mutated_seq.push_back(new_seq[i]);
        }
    }
    return char_edits + max_blocks + reverses;
}

/*
./rearrangement_sim filename no_samples m_min m_max l_min l_max move_remove_rate max_iterations char_error_rate step_interval

Example call:
./rearrangement_sim viral.fna 10 100 150 15 30 3 5 10 3
*/
int main(int argc, char *argv[]) { 

    if (argc != 11) {
        printf("Wrong number of arguments.\n");
        printf("Please run the rearrangement simulator by ./rearrangement_sim filename no_samples m_min m_max l_min l_max move_remove_rate max_iterations error_rate step_interval\n");
        return 0;
    }

    auto start = std::chrono::high_resolution_clock::now();

    char *filename = argv[1];    

    int samples = atoi(argv[2]);
    int m_min = atoi(argv[3]);
    int m_max = atoi(argv[4]);
    int l_min = atoi(argv[5]);
    int l_max = atoi(argv[6]);
    int move_remove_rate = atoi(argv[7]);
    int max_iterations = atoi(argv[8]);

    double char_error_rate = atoi(argv[9]) / 100.0;
    int intensity_test = atoi(argv[10]);

    double reverse_rate = 0.15;

    int m_range = m_max - m_min + 1;
    
    printf("Running %d simulations on sequence length range: (%d, %d); and block length range: (%d, %d)\n", samples * (90 / intensity_test), m_min, m_max, l_min, l_max);

    std::string s1;
    std::string s;
    std::ifstream f;
    f.open(filename);
    std::getline(f, s1);
    for (int i = 0; i < 350; i++) {
        std::getline(f, s1);
        s.append(s1);
    }
    f.close();
    int len = s.size();

    int total_range = 0;
    int total_miss = 0;

    timeval t1, t2;
    double elapsed;


    for (int j = 10; j <= 100; j += intensity_test) {
        double intensity = j / 100.0;
        printf("intensity: %.2f\n", intensity);
        for (int i = 0; i < samples; i++) {
            int m = m_min + (std::rand() % m_range);
            int start = std::rand() % (len - m);
            
            std::string seq = s.substr(start, m);
            std::string mutated_seq;


            int max_block;
            
            max_block = std::floor(((m + 0.0) / (l_max + 0.0)) * intensity);

            int removes = max_block / (move_remove_rate+1);
            int moves = max_block - removes;

            int simulated_bed = mutate_seq(mutated_seq, seq, reverse_rate, char_error_rate, l_min, l_max, removes, moves);
            std::vector<struct block_match> block_matches;
            int n = mutated_seq.size();
            printf("m = %d, n = %d.", m, n);
            fflush(stdout);

            auto t1 = std::chrono::high_resolution_clock::now();
            int calculated_bed = bed(block_matches, seq, mutated_seq, max_iterations, l_min, l_max, 2 * char_error_rate);
            int calculated_blocks = block_matches.size(); 
            auto t2 = std::chrono::high_resolution_clock::now();
            
            int pairwise = lev_dist_edlib(seq, mutated_seq, INF_DIST);
            double miss_rate = (calculated_bed - simulated_bed + 0.0) / (pairwise - simulated_bed + 0.0);

            std::chrono::duration<double> elapsed = t2 - t1;

            printf("\tTook %.2f seconds\n", elapsed.count());
            printf("\tSimulated: %d, calculated: %d, pairwise: %d. miss-rate = %.2f%%.\n", simulated_bed, calculated_bed, pairwise, miss_rate * 100);
            printf("\tSimulated block count: %d, calculated block count: %d\n\n", max_block, calculated_blocks);

            total_miss += calculated_bed - simulated_bed;
            total_range += pairwise - simulated_bed;
        }

        double avg_miss_rate = (total_miss + 0.0) / (total_range + 0.0);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end - start;

        printf("For intensity %d, over %d simulations, overall miss rate = %.3f%%, and accuracy = %.3f%%.\n", j, samples, avg_miss_rate * 100, (1-avg_miss_rate) * 100);
        printf("\ttotal_miss = %d, total_range = %d\n", total_miss, total_range);
        printf("Simulations took %.1f seconds\n\n", elapsed.count());
    }
}
