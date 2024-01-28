#include "bed.h"   
#include "utils.h"
#include "constants.h"

#include <chrono>
#include <getopt.h>
#include <iostream>
#include <string>
#include <cstring>


#define DEFAULT_L_MIN 8
#define DEFAULT_L_MAX 15
#define DEFAULT_MAX_ITERATIONS 3
#define DEFAULT_ERROR_RATE 0.3

#define VERSION "1.2"
#define MAXLINE 32

/*
void print_W(std::vector<std::vector<int>> &W) {
    for (std::size_t i = 0; i < W.size(); i++) {
        for (std::size_t j = 0; j < W[i].size(); j++) {
            if (W[i][j] != INF_DIST) {
                std::cout << "W[" << i << ", " << j << "] = " << W[i][j] << std::endl;
            }
        }
    }
    std::cout << std::endl;
}

void print_S(std::vector<std::vector<std::tuple<int, int, bool>> > &S) {
    for (std::size_t i = 0; i < S.size(); i++) {
        for (std::size_t j = 0; j < S[i].size(); j++) {
            if (std::get<1>(S[i][j]) != -1) {
                std::cout << "S[" << i << ", " << j << "] =  (" << std::get<0>(S[i][j]) << ", " << std::get<1>(S[i][j]) 
                    << ", " << std::get<2>(S[i][j]) << ")" << std::endl;
            }
        }
    }
    std::cout << std::endl;
}

void print_N(std::vector<int> &N) {
    for (std::size_t i = 0; i < N.size(); i++) {
            std::cout << "N[" << i << "] = " << N[i] << std::endl;
    }
}
*/
void print_block_matches(std::vector<struct block_match> &matches, dna_sequence &seq1, dna_sequence &seq2, FILE *out) {
    fprintf(out, "\n");
    for (std::size_t i = 0; i < matches.size(); i++) {
        int i1, j1, i2, j2, reversed, distance;
        i1 = matches[i].loc1.first;
        j1 = matches[i].loc1.second;
        i2 = matches[i].loc2.first;
        j2 = matches[i].loc2.second;
        reversed = matches[i].reversed; 
        distance = matches[i].dist;

        dna_block b1(&seq1, i1, j1);
        
        if (matches[i].remove) {
            fprintf(out, "Block remove S[%d, %d)\n", i1, j1);            
            fprintf(out, "\t%s\n", b1.c_str());
        } else {
            dna_block b2(&seq2, i2, j2);        

            fprintf(out, "Block move S[%d, %d) to T[%d, %d). distance: %d", i1, j1, i2, j2, distance - 1);
            if (reversed)
                fprintf(out, " (reversed)");
            
            fprintf(out, "\n");
            fprintf(out, "\t%s\n", b1.c_str());
            fprintf(out, "\t%s\n", b2.c_str());
        }
        fprintf(out, "\n"); 
    }
}

void help() {
    printf("Usage: \n");
    printf("./saber -s source.fa -t target.fa [-optional arguments]\n\n");

    printf("Arguments:\n\n");

    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "Name","(-code)", "Required", "Type", "Description");
    printf("--------------------------------------------------------------------------------\n");
    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "version", "(-v)", "no", "Flag", "Display the current version of SABER.");
    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "help", "(-h)", "no", "Flag", "Display help screen");
    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "source", "(-s)", "yes", "String", "Specify the file containing the source sequence (must be fasta or fastq file)");
    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "target", "(-t)", "yes", "String", "Specify the file containing the target sequence (must be fasta or fastq file)");
    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "out", "(-o)", "no", "String", "Specify the output filename (Default: stdout)");
    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "iterations", "(-i)", "no", "Integer", "Specify the maximum number of iterations (Default: 1)");
    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "l-min", "(-l)", "no", "Integer", "Specify the minimum block length (Default: 8)");
    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "l-max", "(-m)", "no", "Integer", "Specify the maximum block length (Default: 15)");
    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "error", "(-e)", "no", "Float", "Specify the error rate (Default: 0.3)");
    printf("%-12s  %-9s  %-9s  %-8s  %s\n", "runtime", "(-r)", "no", "Flag", "Display time taken during the computation");
}

void error(const char * msg) {
    fprintf(stderr, "%s\n", msg);
    fprintf(stderr, "./saber --help for more information\n");
    exit(0);
}

void version() {
    printf("SABER %s\n", VERSION);
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
       return dot + 1;
}

bool is_fasta(char *filename) {
    
    char *ext = strrchr(filename, '.');
    if(!ext || ext == filename) return false;
    
    if (strcasecmp(ext, "fa") || strcasecmp(ext, "fasta") || strcasecmp(ext, "fq") || strcasecmp(ext, "fastq"))
        return true;

    return false;
}


// TODO: add dont-report-sequences argument 
int main(int argc, char ** argv) {
    srand((time(0)));

    // Read the arguments to set parameter values
    // help, source, target, output, max_iterations, l-min, l-max, errorrate 

    char s[MAXLINE] = "";
    char t[MAXLINE] = "";
    char o[MAXLINE] = "";

    FILE *out = stdout;

    bool output_set = false;
    bool report_time = false;
    
    int max_iterations = DEFAULT_MAX_ITERATIONS;
    int min_block = DEFAULT_L_MIN; 
    int max_block = DEFAULT_L_MAX;
    double error_rate = DEFAULT_ERROR_RATE;

    bool i_set = false, l_set = false, m_set = false, e_set = false;


    struct option long_opt[] = {
        {"version", 0, NULL, 'v'},
        {"help", 0, NULL, 'h'},
        {"out", 0, NULL, 'o'},
        {"source", 1, NULL, 's'},
        {"target", 1, NULL, 't'},
        {"iterations", 0, NULL, 'i'},
        {"l_min", 0, NULL, 'l'},
        {"l_max", 0, NULL, 'm'},
        {"error", 0, NULL, 'e'},
        {"runtime", 0, NULL, 'r'},
        {NULL, 0, NULL, 0}
    };

    char opt;
    do {
        opt = getopt_long(argc, argv, "vho:s:t:i:l:m:e:r", long_opt, NULL);

        switch (opt) {
        case 'v':
            version();
            exit(0);
        case 'h':
            help();
            exit(0);
        case 'o':
            strncpy(o, optarg, MAXLINE);
            out = fopen(o, "w");
            output_set = true;
            break;
        case 's':
            strncpy(s, optarg, MAXLINE);
            break;
        case 't':
            strncpy(t, optarg, MAXLINE);
            break;
        case 'i':
            max_iterations = atoi(optarg);
            i_set = true;
            break;
        case 'l':
            min_block = atoi(optarg);
            l_set = true;
            break;
        case 'm':
            max_block = atoi(optarg);
            m_set = true;
            break;
        case 'e':
            
            error_rate = atof(optarg);
            e_set = true;
            break;
        case 'r':
            report_time = true;
            break;
        }

    } while (opt != -1);

    // Check the validity of arguments
    if (strlen(s) == 0) 
        error("Source filename must be specified.");
    
    if (strlen(t) == 0) 
        error("Target filename must be specified.");
    
    if (!is_fasta(s) || !is_fasta(t)) 
        error("Only fasta and fastq files are supported.");
    
    if (max_block < min_block)
        error("Minimum block length cannot be larger than maximum block length");

    if (error_rate > 1)
        error("Error rate must not exceed 1.");

    if (error_rate < 0)
        error("Error rate cannot be negative.");
    
    if (max_iterations <= 0) 
        error("The number of iterations must be greater than or equal to 1");
    
    if (max_block <= 0 || min_block <= 0)
        error("Block length limits must be positive integers.");

    // More checks
    if (!output_set)
        printf("Output path was not specified, printing the output directly to the console...\n");

    
    if (!l_set && !m_set) {
        printf("Block lengths were not specified, using the default values...\n");

    } else if (!m_set) {
        printf("Minimum block length was specified but maximum block length was not specified. Maximum block length was automatically set...\n");
        max_block = min_block * 2 - 1;

    }  else if (!l_set) {
        printf("Maximum block length was specified but minimum block length was not specified. Minimum block length was automatically set...\n");
        min_block = (max_block + 1) / 2;

    }

    if (!e_set) 
        printf("Error-rate was not specified, using the default values...\n");
    
    if (!i_set) 
        printf("Number of iterations was not specified, using the default values...\n");

    fprintf(out, "l-min = %d, l-max = %d, max-iterations = %d, error-rate = %.2f\n", min_block, max_block, max_iterations, error_rate);

    // Read fasta or fastq files to get two sequences
    dna_sequence seq1, seq2;

    if (!read_sequence(s, seq1)) {
        std::cerr << "Failed to read the first sequence from file " << s << std::endl;
        return 1;
    }

    if (!read_sequence(t, seq2)) {
        std::cerr << "Failed to read the first sequence from file " << t << std::endl;
        return 1;
    }
 
    std::vector<std::vector<int>> W;
    std::vector<std::vector<std::tuple<int,int,bool>> > S;
    std::vector<int> N;
    std::vector<struct block_match> block_matches;
    std::string new_seq1, new_seq2;

    std::chrono::time_point<std::chrono::high_resolution_clock> t1, t2;
    double total_elapsed = 0;

    if (report_time)
        t1 = std::chrono::high_resolution_clock::now();

    calculate_W(W, S, seq1, seq2, min_block, max_block, error_rate);    
    if (report_time) {
        t2 = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = t2 - t1;
        total_elapsed += elapsed.count();

        fprintf(out, "\nW computed in %.2f seconds.\n", elapsed.count());

        t1 = std::chrono::high_resolution_clock::now();
    }

    
    calculate_N(N, W, S, seq1, seq2, max_iterations, min_block, max_block);

    //print_N(N);
    
    if (report_time) {
        t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = t2 - t1;
        total_elapsed += elapsed.count();
        fprintf(out, "N computed in %.2f seconds.\n", elapsed.count());
    }

    fprintf(out, "\n");
    int block_edit_distance = block_edit_score(block_matches, new_seq1, new_seq2, seq1, seq2, N, N.size() - 1, W, S, INF_DIST, min_block, max_block);

    fprintf(out, "Block Operations: \n");
    fprintf(out, "-------------------------------------\n");
    print_block_matches(block_matches, seq1, seq2, out);
    fprintf(out, "-------------------------------------\n");

    fprintf(out, "\nRemaining sequences:\n\n");
    
    fprintf(out, "%s\n", new_seq1.c_str());
    fprintf(out, "%s\n", new_seq2.c_str());
    fprintf(out, "\nDistance of remaining sequences: %d\n\n", lev_dist_edlib(new_seq1, new_seq2, INF_DIST));
    
    fprintf(out, "\nBED(S, T) = %d\n", block_edit_distance);

    if (report_time) {
        fprintf(out, "\nThe algorithm took %.2f seconds\n", total_elapsed);
    }

    fclose(out);

    return 0;
} 