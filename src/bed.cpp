#include "bed.h"
#include "constants.h"
#include "utils.h"

void construct_match_table(std::vector<std::vector<int>> &d, dna_block &_block, dna_sequence &_B) {
	unsigned int l = _block.size();
	unsigned int n = _B.size();

	d = std::vector<std::vector<int>>(l + 1, std::vector<int>(n + 1, 0));

	// initial conditions
	for(unsigned int i = 1; i <= l; ++i) d[i][0] = i * C_CHAR_DELETION;

	// computation of the recursive formulation
	for(unsigned int i = 1; i <= l; ++i)
		for(unsigned int j = 1; j <= n; ++j)
            d[i][j] = std::min({d[i - 1][j] + C_CHAR_DELETION, d[i][j - 1] + C_CHAR_INSERTION, d[i - 1][j - 1] + (_block[i - 1] == _B[j - 1] ? 0 : C_CHAR_SNP)});
}


int find_starting_position(std::vector<std::vector<int>> &d, int _j, int best_end, dna_block &_block, dna_sequence &_B) {
	int j = best_end;
	int i = _j;

	while (i > 0 && j > 0) {
		if (d[i-1][j] + C_CHAR_DELETION == d[i][j]) { // Deletion
			i--;
		} else if (_block[i-1] == _B[j-1] && d[i-1][j-1] == d[i][j]) { // Character match
			i--;
			j--;
		} else if (_block[i-1] != _B[j-1] && d[i-1][j-1] + C_CHAR_SNP == d[i][j]) { // SNP
			i--;
			j--;
        } else if (d[i][j-1] + C_CHAR_INSERTION == d[i][j]) { // Insertions
			j--;
		} else {
		    throw std::runtime_error("Invalid match table.");
		}
	}

	return j;
}


bool ignore_match(int dist, int l1, int l2, double error_rate) {
    return dist > std::ceil(error_rate * ((l1 + l2 + 0.0) / 2.0));
}

void calculate_W(std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, dna_sequence &seq1, dna_sequence &seq2, int min_block, int max_block, double error_rate) {
    // Random initialization
    coin_tosser ct;

    int m = seq1.size();
    int n = seq2.size();

    int block_range = max_block - min_block + 1;

    W = std::vector<std::vector<int>>(m, std::vector<int>(block_range, INF_DIST));
    S = std::vector<std::vector<std::tuple<int,int, bool>> >(m, std::vector<std::tuple<int,int,bool>>(block_range, UNKNOWN_TUPLE));

    for (int i = 0; i < m; i++) {
        int l_max;

	    if (max_block + i < m)
		    l_max = max_block;
	    else
		    l_max = m - i;
	    

        if (l_max < min_block) continue;

        std::vector<std::vector<int>> d;
	    std::vector<std::vector<int>> d_rev;

	    dna_block block = dna_block(&seq1, i, i + l_max);
        dna_block block_rev = dna_block(&seq1, i, i + l_max, true);

	    construct_match_table(d, block, seq2);
	    construct_match_table(d_rev, block_rev, seq2);

	    for (int j = l_max; j >= min_block; --j) {
		    int j_idx = j - min_block;
		

		    // To find the best match for block of length j using the match table,
		    // 		Find min from both tables for row j.
		    // 		Determine whether to reverse or not by selecting the best table for this j.
		    // 		Find the starting point of the block by going backwards from the table.

		    int dist = INF_DIST;
		    int rev_dist = INF_DIST;

		    int loc = -1;
		    int rev_loc = -1; 

		    for (int k = 0; k <= n; k++) {
			    if (d[j][k] < dist || (d[j][k] == dist && ct.toss())) {
				    dist = d[j][k];
				    loc = k;
			    } 

			    if (d_rev[j][k] < rev_dist || (d_rev[j][k] == rev_dist && ct.toss())) {
				    rev_dist = d_rev[j][k];
				    rev_loc = k;
			    }
 		    }

		    bool reversed = rev_dist < dist;
		    int best_dist, best_end, best_start;

		    if (reversed) {
			    // Find starting position from d_rev
			    best_end = rev_loc;
			    best_dist = rev_dist + C_BLOCK_MOVE + C_BLOCK_REVERSE;
			    best_start = find_starting_position(d_rev, j, best_end, block_rev, seq2);

		    } else {
			    // Find starting position from d 
			    best_end = loc;
			    best_dist = dist + C_BLOCK_MOVE;
			    best_start = find_starting_position(d, j, best_end, block, seq2);
		    }

            if (best_start == -1 || ignore_match(best_dist, j, best_end - best_start, error_rate)) continue;
            
			S[i][j_idx] = std::make_tuple(best_start, best_end, reversed);
            W[i][j_idx] = best_dist;
	    }
    }
}

void calculate_N(std::vector<int> &N, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, dna_sequence &seq1, dna_sequence &seq2, int max_iterations, int min_block, int max_block) {
    int m = seq1.size();

    N = std::vector<int>(m+1, -1);

    int char_dist = lev_dist_edlib(seq1, seq2, INF_DIST);

    int bs = char_dist;

    std::vector<int> order(m - min_block + 1);
    std::iota(order.begin(), order.end(), min_block);

    for (int it = 0; it < max_iterations; it++) {
        // For the first iteration, calculate normally;
        // In the next iterations, calculate according to a randomized order.

        if (it > 0) {
            std::random_shuffle(order.begin(), order.end());
        }

        int score_before = bs;
        for (int i : order) { 
            int choice_before = N[i];
            
            int best_choice = N[i];
            int best_score = bs;
            
            // pass: dont consider this character in a block operation
            N[i] = -1;
            int score = block_edit_score(seq1, seq2, N, m, W, S, best_score - 1, min_block, max_block);
            if (score <= best_score - 1) {
                best_score = score;
                best_choice = -1;
            }

            int dont_remove = 0;

            for (int j = max_block; j >= min_block; j--) {
                
                int i_idx = i-j+1;
                int j_idx = j-min_block;
                
                if (i_idx < 0)
                    continue;

                // Block Move
                if (W[i_idx][j_idx] != INF_DIST && !block_overlaps_in_B(N, m, S, min_block, max_block)) {

                    N[i] = j;
                    score =  block_edit_score(seq1, seq2, N, m, W, S, best_score - 1, min_block, max_block);
                    if (score < best_score) {
                        best_score = score;
                        best_choice = j;
                    }
                }
                

                // Block remove
                if (dont_remove <= 0) {

                    N[i] = -j;
                    score =  block_edit_score(seq1, seq2, N, m, W, S, INF_DIST, min_block, max_block);
                    if (score < best_score) {
                        best_score = score;
                        best_choice = -j;
                    } else {
                        dont_remove = score - best_score - 1;
                    }
                }
                dont_remove--;
            }
            
            if (best_score < bs) {
                N[i] = best_choice;
                bs = best_score;
            } else {
                N[i] = choice_before;
            }
        }

        if (bs == score_before) break;
    }
}

void get_block_matches(std::vector<struct block_match> &block_matches, std::vector<int> &N, int i, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int min_block, int max_block) {
    while (i >= min_block) {
        if (N[i] == -1) {
            // Pass Character
            i--;

        } else if (N[i] < -1) {
            // Block Remove
            int j = -N[i];
            
            if (i-j+1 < 0)
                break;

            struct block_match match = {
                .loc1 = std::make_pair(i-j+1, i+1),
                .remove = true
            };

            block_matches.emplace_back(match);         

            i -= j;
        } else {
            // Block Move or Reverse
            int j = N[i];
            
            if (i-j+1 < 0)
                break;

            struct block_match match = {
                .loc1 = std::make_pair(i-j+1, i+1),
                .loc2 = std::make_pair(std::get<0>(S[i-j+1][j-min_block]), std::get<1>(S[i-j+1][j-min_block])),
                .remove = false,
                .reversed = std::get<2>(S[i-j+1][j-min_block]),
                .dist = W[i-j+1][j-min_block]
            };

            block_matches.emplace_back(match);         

            i -= j; 
        }
    }
}

void get_remaining_characters(std::string &new_seq1, std::string &new_seq2, dna_sequence &seq1, dna_sequence &seq2, int _i, std::vector<struct block_match> &block_matches) {
    int n = seq2.size();
    
    new_seq1 = "";
    new_seq2 = "";

    for (int i = 0; i < _i; i++) {
        bool occupied = false;
        for (struct block_match match : block_matches) {
            if (i >= match.loc1.first && i < match.loc1.second) {
                occupied = true;
                break;
            }
        }
        if (!occupied) {
            new_seq1.push_back(seq1[i]);
            //std::cout << "seq1: unoccupied char at " << i << ": " << seq1[i] << std::endl;
        }
    }

    for (int i = 0; i < n; i++) {
        bool occupied = false;
        for (struct block_match match : block_matches) {
            if (!match.remove && i >= match.loc2.first && i < match.loc2.second) {
                occupied = true;
                break;
            }
        }
        if (!occupied) {
            new_seq2.push_back(seq2[i]);
        }
    }
}

void get_remaining_characters(std::string &new_seq1, std::string &new_seq2, dna_sequence &seq1, dna_sequence &seq2, std::vector<bool> &marked1, std::vector<bool> &marked2) {
    int m = seq1.size();
    int n = seq2.size();
    
    new_seq1 = "";
    new_seq2 = "";

    for (int i = 0; i < m; i++) 
        if (!marked1[i])
            new_seq1.push_back(seq1[i]);
    
    for (int j = 0; j < n; j++) 
        if (!marked2[j])
            new_seq2.push_back(seq2[j]);
}

bool overlap(int b1_s, int b1_e, int b2_s, int b2_e) {
    bool dont_overlap = b1_s >= b2_e || b2_s >= b1_e;
    return !dont_overlap;
}

bool block_overlaps_in_B(std::vector<int> &N, int _i, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int min_block, int max_block) {
    // precondition: N[i] >= MIN_BLOCK_LEN
    if (N[_i] < min_block) 
        return false;
    
    int i = _i;
    int j = N[i]; 
    int i_idx = i - j + 1;
    int j_idx = j - min_block;
    
    std::tuple<int,int, bool> block = S[i_idx][j_idx];
    int b_start = std::get<0>(block);
    int b_end = std::get<1>(block);

    i -= j;
    while (i >= min_block - 1) {
        if (N[i] < 0) { // character pass or block remove, overlap is not possible
            j = -N[i];
            i -= j;
        } else {
            j = N[i];
            i_idx = i - j + 1;
            j_idx = j - min_block;

            std::tuple<int,int, bool> block2 = S[i_idx][j_idx];
            int b2_start = std::get<0>(block2);
            int b2_end = std::get<1>(block2);

            if (overlap(b_start, b_end, b2_start, b2_end))
                return true;
            i -= j;       
        }   
    }
    return false;
}

// Returns the block edit score as well as block matches and the new sequences
int block_edit_score(std::vector<struct block_match> &block_matches, std::string &new_seq1, std::string &new_seq2, dna_sequence &seq1, dna_sequence &seq2, std::vector<int> &N, int _i, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int cutoff, int min_block, int max_block) {
    block_matches = std::vector<struct block_match>();
    int m = seq1.size();
    int n = seq2.size();

    get_block_matches(block_matches, N, _i, W, S, min_block, max_block);
    
    // Calculate the block score and the remaining characters at the same time
    std::vector<bool> marked1(m, false), marked2(n, false);

    int block_dist = 0;
    for (struct block_match match : block_matches) {
        for (int i = match.loc1.first; i < match.loc1.second; i++)   marked1[i] = true;

        if (match.remove) {
            block_dist += C_BLOCK_DELETION;
            
        } else {
            block_dist += match.dist;
            for (int j = match.loc2.first; j < match.loc2.second; j++)   marked2[j] = true;
        }
    }

    //std::cout << "\t\tblock dist = " << block_dist << std::endl;

    int char_cutoff = cutoff - block_dist;
    if (char_cutoff < 0) return INF_DIST;
    //std::cout << "\t\tchar cutoff = " << char_cutoff << std::endl;
    get_remaining_characters(new_seq1, new_seq2, seq1, seq2, marked1, marked2);

    int remaining_char_dist = lev_dist_edlib(new_seq1, new_seq2, char_cutoff);
    //std::cout << "\t\tchar dist = " << remaining_char_dist << std::endl;

    if (remaining_char_dist == -1 || remaining_char_dist > char_cutoff) return INF_DIST;

    return block_dist + remaining_char_dist;
}

// Does not return block matches or new sequences
int block_edit_score(dna_sequence &seq1, dna_sequence &seq2, std::vector<int> &N, int i, std::vector<std::vector<int>> &W, std::vector<std::vector<std::tuple<int,int,bool>> > &S, int cutoff, int min_block, int max_block) {
    std::vector<struct block_match> block_matches;
    std::string new_seq1, new_seq2;
    return block_edit_score(block_matches, new_seq1, new_seq2, seq1, seq2, N, i, W, S, cutoff, min_block, max_block);
}

// Precondition: block matches were already computed
int block_edit_score(dna_sequence &seq1, dna_sequence &seq2, std::vector<struct block_match> &block_matches) {
    int m = seq1.size();

    std::string new_seq1, new_seq2;
    get_remaining_characters(new_seq1, new_seq2, seq1, seq2, m, block_matches);

    int block_dist = 0;
    for (struct block_match match : block_matches) {
        if (match.remove) {
            block_dist += C_BLOCK_DELETION;
        } else {
            block_dist += match.dist;
        }
    }

    int remaining_char_dist = lev_dist_edlib(new_seq1, new_seq2, -1);

    return block_dist + remaining_char_dist;
}


int bed(std::vector<struct block_match> &block_matches, std::string &_seq1, std::string &_seq2, int max_iterations, int min_block, int max_block, double error_rate) {
    std::vector<std::vector<int>> W;
    std::vector<std::vector<std::tuple<int,int,bool>> > S;
    std::vector<int> M;
    std::vector<int> N;

    std::string new_seq1, new_seq2;

    calculate_W(W, S, _seq1, _seq2, min_block, max_block, error_rate);
    calculate_N(N, W, S, _seq1, _seq2, max_iterations, min_block, max_block);

    int block_edit_distance = block_edit_score(block_matches, new_seq1, new_seq2, _seq1, _seq2, N, N.size() - 1, W, S, INF_DIST, min_block, max_block);
    return block_edit_distance;
}
