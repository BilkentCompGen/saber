#include "edit_dist.h"

// Standard Levenshtein Distance
int lev_dist(std::string &_a, std::string &_b) {
    int m = _a.length(), n = _b.length();
    int d[2][n + 1];
    d[0][0] = 0;
    d[1][0] = 1;
    for (int i = 0; i < n + 1; i++) d[0][i] = i;
    for (int i = 1; i < m + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            d[i&1][j] = min(min(d[(i-1)&1][j], d[i&1][j-1]) + 1, d[(i-1)&1][j-1] + (_a[i-1] == _b[j-1] ? 0 : 1));
        }
    }
    return d[m&1][n];
}

int lev_dist(dna_block &_a, dna_block &_b) {
    int m = seqan::length(_a), n = seqan::length(_b);
    int d[2][n + 1];
    d[0][0] = 0;
    d[1][0] = 1;
    for (int i = 0; i < n + 1; i++) d[0][i] = i;
    for (int i = 1; i < m + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            d[i&1][j] = min(min(d[(i-1)&1][j], d[i&1][j-1]) + 1, d[(i-1)&1][j-1] + (_a[i-1] == _b[j-1] ? 0 : 1));
        }
    }
    return d[m&1][n];
}

int lev_dist(dna_block &_a, std::string &_b) {
    int m = seqan::length(_a), n = _b.length();
    int d[2][n + 1];
    d[0][0] = 0;
    d[1][0] = 1;
    for (int i = 0; i < n + 1; i++) d[0][i] = i;
    for (int i = 1; i < m + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            d[i&1][j] = min(min(d[(i-1)&1][j], d[i&1][j-1]) + 1, d[(i-1)&1][j-1] + (_a[i-1] == _b[j-1] ? 0 : 1));
        }
    }
    return d[m&1][n];
}

int lev_dist(dna_block &_a, seqan::DnaString &_b) {
    int m = seqan::length(_a), n = seqan::length(_b);
    uint32_t d[2][n + 1];
    d[0][0] = 0;
    d[1][0] = 1;
    for (int i = 0; i < n + 1; i++) d[0][i] = i;
    for (int i = 1; i < m + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            d[i&1][j] = min(min(d[(i-1)&1][j], d[i&1][j-1]) + 1, d[(i-1)&1][j-1] + (_a[i-1] == _b[j-1] ? 0 : 1));
        }
    }
    return d[m&1][n];
}

int lev_dist(seqan::DnaString &_a, seqan::DnaString &_b) {
    int m = seqan::length(_a), n = seqan::length(_b);
    uint32_t d[2][n + 1];
    d[0][0] = 0;
    d[1][0] = 1;
    for (int i = 0; i < n + 1; i++) d[0][i] = i;
    for (int i = 1; i < m + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            d[i&1][j] = min(min(d[(i-1)&1][j], d[i&1][j-1]) + 1, d[(i-1)&1][j-1] + (_a[i-1] == _b[j-1] ? 0 : 1));
        }
    }
    return d[m&1][n];
}


int lev_dist_edlib(seqan::DnaString &_a, seqan::DnaString &_b, int k) {
	std::string a = to_string(_a);
	std::string b = to_string(_b);

	return lev_dist_edlib(a, b, k);
}

int lev_dist_edlib(std::string &_a, std::string &_b, int k) {
	int m = _a.length(), n = _b.length();

	EdlibEqualityPair *eep = NULL;
	EdlibAlignConfig config = edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, eep, 0);

	EdlibAlignResult result = edlibAlign(_a.c_str(), m, _b.c_str(), n, config);  	

	int dist = result.editDistance;
	edlibFreeAlignResult(result);
	return dist;
}

int i_base(char c) {
	switch (c) 
	{
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	default:
		return 0;
	}
}

int lev_dist_myers_64(std::string &_a, std::string &_b) {
	int m = _a.size(), n = _b.size();

    uint64_t Peq[4];
    uint64_t Eq, Xv, Xh, Ph, Mh, Pv, Mv, Last;
    int64_t i;
    int Score = n;

    memset(Peq, 0, sizeof(Peq));

    for (i = 0; i < n; i++)
        Peq[i_base(_b[i])] |= (uint64_t) 1 << i;

    Mv = 0;
    Pv = (uint64_t) -1;
    Last = (uint64_t) 1 << (n - 1);

    for (i = 0; i < m; i++) {
        Eq = Peq[i_base(_a[i])];

        Xv = Eq | Mv;
        Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

        Ph = Mv | ~ (Xh | Pv);
        Mh = Pv & Xh;

        if (Ph & Last) Score++;
        if (Mh & Last) Score--;

        Ph = (Ph << 1) | 1;
        Mh = (Mh << 1);

        Pv = Mh | ~ (Xv | Ph);
        Mv = Ph & Xv;
    }
    return Score;
}


std::string to_string(seqan::DnaString &seq) {
	int n = seqan::length(seq);

	std::string result;

	for(int i = 0; i < n; ++i) 
		result.push_back(seq[i]);	
	
	return result;
}

std::string to_string(dna_block &seq) {
	int n = seqan::length(seq);

	std::string result;

	for(int i = 0; i < n; ++i) 
		result.push_back(seq[i]);	
	
	return result;
}

// Reverse Complement Operations
std::string reverse_complement(seqan::DnaString &seq) {
	int n = seqan::length(seq);

	std::string result;

	for (int i = n-1; i >= 0; i--) {
		if (seq[i] == 'A') 
			result.push_back('T');
		else if (seq[i] == 'T')
			result.push_back('A');
		else if (seq[i] == 'C')
			result.push_back('G');
		else if (seq[i] == 'G')
			result.push_back('C');
	}

	return result;
}

std::string reverse_complement(dna_block &seq) {
	int n = seqan::length(seq);

	std::string result;

	for (int i = n-1; i >= 0; i--) {
		if (seq[i] == 'A') 
			result.push_back('T');
		else if (seq[i] == 'T')
			result.push_back('A');
		else if (seq[i] == 'C')
			result.push_back('G');
		else if (seq[i] == 'G')
			result.push_back('C');
	}

	return result;
}

std::string reverse_complement(std::string &seq) {
	int n = seq.size();

	std::string result;

	for (int i = n-1; i >= 0; i--) {
		if (seq[i] == 'A') 
			result.push_back('T');
		else if (seq[i] == 'T')
			result.push_back('A');
		else if (seq[i] == 'C')
			result.push_back('G');
		else if (seq[i] == 'G')
			result.push_back('C');
	}
	
	return result;
}

int dist_rev(bool &_reversed, std::string &_a, std::string &_b) {
	std::string b_rev = reverse_complement(_b);
	int distance = lev_dist(_a, _b);
	int distance_rev = lev_dist(_a, b_rev) + C_BLOCK_REVERSE;
	if (distance <= distance_rev) {
		_reversed = false;
		return distance;
	} else {
		_reversed = true;
		return distance_rev;
	}
}

int dist_rev(bool &_reversed, dna_block &_a, dna_block &_b, std::string &_b_rev) {
	
	int distance = lev_dist(_a, _b);
	int distance_rev = lev_dist(_a, _b_rev) + C_BLOCK_REVERSE;

	if (distance <= distance_rev) {
		_reversed = false;
		return distance;
	} else {
		_reversed = true;
		return distance_rev;
	}
}

int dist_rev(bool &_reversed, dna_block &_a, dna_block &_b) {
	std::string b_rev = reverse_complement(_b);

	int distance = lev_dist(_a, _b);
	int distance_rev = lev_dist(_a, b_rev) + C_BLOCK_REVERSE;

	if (distance <= distance_rev) {
		_reversed = false;
		return distance;
	} else {
		_reversed = true;
		return distance_rev;
	}

}

int dist_rev(bool &_reversed, seqan::DnaString &_a, seqan::DnaString &_b) {
	seqan::DnaString b_rev = reverse_complement(_b);
	int distance = lev_dist(_a, _b);
	int distance_rev = lev_dist(_a, b_rev) + C_BLOCK_REVERSE;
	if (distance <= distance_rev) {
		_reversed = false;
		return distance;
	} else {
		_reversed = true;
		return distance_rev;
	}
}

