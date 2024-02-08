#include "bed.h"
#include "constants.h"
#include "utils.h"

#include <cmath>
#include <stack>

void construct_match_table(std::vector<std::vector<int>> &d, dna_block &_block, dna_sequence &_B) {
	unsigned int l = _block.size();
	unsigned int n = _B.size();

	d = std::vector<std::vector<int>>(l + 1, std::vector<int>(n + 1, 0));

	// initial conditions
	for(unsigned int i = 1; i <= l; ++i) d[i][0] = i * C_CHAR_DELETION;

	// computation of the recursive formulation
	for(unsigned int i = 1; i <= l; ++i)
		for(unsigned int j = 1; j <= n; ++j)
            d[i][j] = min_3(d[i-1][j] + C_CHAR_DELETION, d[i][j-1] + C_CHAR_INSERTION, d[i-1][j-1] + (_block[i-1] == _B[j-1] ? 0 : C_CHAR_SNP));
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
    return dist > (error_rate * ((l1 + l2) / 2.0));
}

void calculate_W(std::vector<std::vector<w_entry>> &W, dna_sequence &seq1, dna_sequence &seq2, alignment_setting settings) {
    // Random initialization
    coin_tosser ct;

    int m = seq1.size();
    int n = seq2.size();

    int block_range = settings.max_block - settings.min_block + 1;

    W = std::vector<std::vector<w_entry>>(m, std::vector<w_entry>(block_range, w_entry()));

    for (int i = 0; i < m; i++) {
        int l_max = std::min(settings.max_block, m - i);
	    
        if (l_max < settings.min_block) continue;

        std::vector<std::vector<int>> d;
	    std::vector<std::vector<int>> d_rev;

	    dna_block block = dna_block(&seq1, i, i + l_max);
        dna_block block_rev = dna_block(&seq1, i, i + l_max, true);

	    construct_match_table(d, block, seq2);
	    construct_match_table(d_rev, block_rev, seq2);

	    for (int j = l_max; j >= settings.min_block; --j) {
		    int j_idx = j - settings.min_block;
		

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

            if (best_start == -1 || ignore_match(best_dist, j, best_end - best_start, settings.error_rate)) continue;
            
            W[i][j_idx].set(best_dist, best_start, best_end, reversed);
	    }
    }
}

std::stack<int> reached_list(std::vector<int> &N) {
    int m = N.size() - 1;
    std::stack<int> reached;
    while (m > 0) {
        m -= std::abs(N[m]);
        reached.push(m);
    }
    return reached;
}

std::pair<int, int> optimize_index(std::vector<int> &N, int i, std::vector<std::vector<w_entry>> &W, int old_score, dna_sequence &seq1, dna_sequence &seq2, alignment_setting settings) {
    int choice_before = N[i];
    
    int best_score = old_score;
    int best_choice = choice_before;

    int score;
    
    // pass: dont consider this character in a block operation
    if (choice_before != -1) {
        N[i] = -1;
        score = block_edit_score(seq1, seq2, N, W, best_score - 1, settings);
        if (score < best_score) {
            best_score = score;
            best_choice = -1;
        }
    }

    int dont_remove = 0;

    for (int j = settings.max_block; j >= settings.min_block; j--) {
        int i_idx = i-j+1;
        int j_idx = j-settings.min_block;
                
        if (i_idx < 0) continue;

        // Block Move
        if (W[i_idx][j_idx].valid() && choice_before != j && !block_overlaps_in_B(N, i, j, W, settings)) {
            N[i] = j;
            score =  block_edit_score(seq1, seq2, N, W, best_score - 1, settings);
            if (score < best_score) {
                best_score = score;
                best_choice = j;
            }
        }
                
        // Block remove
        if (dont_remove <= 0 && choice_before != -j) {
            N[i] = -j;
            score =  block_edit_score(seq1, seq2, N, W, INF_DIST, settings);
            if (score < best_score) {
                best_score = score;
                best_choice = -j;
            } else {
                dont_remove = score - best_score - 1;
            }
        }
        dont_remove--;
    }
    N[i] = choice_before;
    
    return std::make_pair(best_choice, best_score);
}

void find_optimal_N(std::vector<int> &N, std::vector<std::vector<w_entry>> &W, dna_sequence &seq1, dna_sequence &seq2, alignment_setting settings) {
    int m = seq1.size();

    N = std::vector<int>(m, -1);

    // Calculate the initial pairwise distance
    int char_dist = lev_dist_edlib(seq1, seq2, INF_DIST);

    // Initial block score
    int bs = char_dist;


    // In the first iteration, go over all indices
    for (int i = settings.min_block; i < m; i++) {
        auto best = optimize_index(N, i, W, bs, seq1, seq2, settings);
        N[i] = best.first;    
        bs = best.second;
    }
    for (int it = 0; it < settings.max_iterations; it++) {
        // Get the list of indices that are reached.
        std::stack<int> reached = reached_list(N);
        
        int score_before = bs;

        while (!reached.empty()) {
            // Get the top index from the stack
            int i = reached.top();
            reached.pop();

            // Optimize the index i and get the best choice and score
            auto best = optimize_index(N, i, W, bs, seq1, seq2, settings);

            // Update N[i] with the best choice
            N[i] = best.first;
            // Update the block score with the best score
            bs = best.second;
        }

        // If the block score didn't improve, break the loop
        if (bs == score_before) break;
    }
}
void get_block_matches(std::vector<block_match> &block_matches, std::vector<int> &N, int i, std::vector<std::vector<w_entry>> &W, alignment_setting settings) {
    while (i >= settings.min_block) {
        if (N[i] == -1) {
            // Pass Character
            i--;

        } else if (N[i] < -1) {
            // Block Remove
            int j = -N[i];
            
            if (i-j+1 < 0)
                break;

            block_match match(i-j+1, i+1);

            block_matches.emplace_back(match);         

            i -= j;
        } else {
            // Block Move or Reverse
            int j = N[i];
            
            if (i-j+1 < 0)
                break;

            block_match match(i-j+1, i+1, W[i-j+1][j-settings.min_block].get_start(), W[i-j+1][j-settings.min_block].get_end(), W[i-j+1][j-settings.min_block].is_reversed(), W[i-j+1][j-settings.min_block].get_score());

            block_matches.emplace_back(match);         

            i -= j; 
        }
    }
}

void get_remaining_characters(dna_sequence &new_seq1, dna_sequence &new_seq2, dna_sequence &seq1, dna_sequence &seq2, std::vector<block_match> &block_matches) {
    int m = seq1.size();
    int n = seq2.size();
    
    new_seq1 = "";
    new_seq2 = "";

    for (int i = 0; i < m; i++) {
        bool occupied = false;
        for (block_match match : block_matches) {
            if (i >= match.loc1.first && i < match.loc1.second) {
                occupied = true;
                break;
            }
        }
        if (!occupied) {
            new_seq1.push_back(seq1[i]);
        }
    }

    for (int i = 0; i < n; i++) {
        bool occupied = false;
        for (block_match match : block_matches) {
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

void get_remaining_characters(dna_sequence &new_seq1, dna_sequence &new_seq2, dna_sequence &seq1, dna_sequence &seq2, std::vector<bool> &marked1, std::vector<bool> &marked2) {
    int m = seq1.size();
    int n = seq2.size();

    new_seq1.reserve(m);
    new_seq2.reserve(n);

    int i = 0;
    for (; i < m && i < n; i++) {
        if (!marked1[i])
            new_seq1.push_back(seq1[i]);
        if (!marked2[i])
            new_seq2.push_back(seq2[i]);
    }

    // If seq1 is longer
    for (; i < m; i++) 
        if (!marked1[i])
            new_seq1.push_back(seq1[i]);
            
    // If seq2 is longer
    for (; i < n; i++) 
        if (!marked2[i])
            new_seq2.push_back(seq2[i]);
}

bool overlap(int b1_s, int b1_e, int b2_s, int b2_e) {
    return (b2_s < b1_e && b2_e > b1_s) && !(b1_s == b2_s && b1_e == b2_e);
}

bool block_overlaps_in_B(std::vector<int> &N, int _i, int _j,std::vector<std::vector<w_entry>> &W, alignment_setting settings) {
    if (_j < settings.min_block || _j > settings.max_block) 
        return false;
    
    int i = _i;
    int j = _j;
    int i_idx = i - j + 1;
    int j_idx = j - settings.min_block;
    
    int b_start = W[i_idx][j_idx].get_start();
    int b_end = W[i_idx][j_idx].get_end();

    std::vector<block_match> block_matches;
    get_block_matches(block_matches, N, N.size() - 1, W, settings);

    for (auto match : block_matches) {
        if (overlap(b_start, b_end, match.loc2.first, match.loc2.second))
            return true;
    }

    return false;
}

// Returns the block edit score as well as block matches and the new sequences
int block_edit_score(std::vector<block_match> &block_matches, dna_sequence &new_seq1, dna_sequence &new_seq2, dna_sequence &seq1, dna_sequence &seq2, std::vector<int> &N, std::vector<std::vector<w_entry>> &W, int cutoff, alignment_setting settings) {
    block_matches = std::vector<block_match>();
    int m = seq1.size();
    int n = seq2.size();

    get_block_matches(block_matches, N, m-1, W, settings);
    
    // Calculate the block score and the remaining characters at the same time
    std::vector<bool> marked1(m, false), marked2(n, false);

    int block_dist = 0;
    for (block_match match : block_matches) {
        for (int i = match.loc1.first; i < match.loc1.second; i++)   marked1[i] = true;

        if (match.remove) {
            block_dist += C_BLOCK_DELETION;
            
        } else {
            block_dist += match.dist;
            for (int j = match.loc2.first; j < match.loc2.second; j++)   marked2[j] = true;
        }

        if (block_dist > cutoff) return INF_DIST;
    }

    int char_cutoff = cutoff - block_dist;
    if (char_cutoff < 0) return INF_DIST;

    get_remaining_characters(new_seq1, new_seq2, seq1, seq2, marked1, marked2);

    int remaining_char_dist = lev_dist_edlib(new_seq1, new_seq2, char_cutoff);

    if (remaining_char_dist == -1 || remaining_char_dist > char_cutoff) return INF_DIST;

    return block_dist + remaining_char_dist;
}

// Does not return block matches or new sequences
int block_edit_score(dna_sequence &seq1, dna_sequence &seq2, std::vector<int> &N, std::vector<std::vector<w_entry>> &W, int cutoff, alignment_setting settings) {
    if (cutoff < 0) return INF_DIST;
    std::vector<block_match> block_matches;
    dna_sequence new_seq1, new_seq2;
    return block_edit_score(block_matches, new_seq1, new_seq2, seq1, seq2, N, W, cutoff, settings);
}

// Precondition: block matches were already computed
int block_edit_score(dna_sequence &seq1, dna_sequence &seq2, std::vector<block_match> &block_matches) {
    dna_sequence new_seq1, new_seq2;
    get_remaining_characters(new_seq1, new_seq2, seq1, seq2, block_matches);

    int block_dist = 0;
    for (block_match match : block_matches) {
        if (match.remove) 
            block_dist += C_BLOCK_DELETION;
        else 
            block_dist += match.dist;
        
    }

    int remaining_char_dist = lev_dist_edlib(new_seq1, new_seq2, -1);

    return block_dist + remaining_char_dist;
}

int bed(std::vector<block_match> &block_matches, dna_sequence &_seq1, dna_sequence &_seq2, int max_iterations, alignment_setting settings) {
    std::vector<std::vector<w_entry>> W;
    std::vector<int> N;

    dna_sequence new_seq1, new_seq2;

    calculate_W(W, _seq1, _seq2, settings);
    find_optimal_N(N, W, _seq1, _seq2, settings);

    int block_edit_distance = block_edit_score(block_matches, new_seq1, new_seq2, _seq1, _seq2, N, W, INF_DIST, settings);
    return block_edit_distance;
}
