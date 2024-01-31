#include "ed.h"
#include "utils.h"

int lev_dist_edlib(dna_block &_a, dna_block &_b, int k) {
	return lev_dist_edlib(_a.c_str(), _a.size(), _b.c_str(), _b.size(), k);
}

int lev_dist_edlib(dna_sequence &_a, dna_sequence &_b, int k) {
	return lev_dist_edlib(_a.c_str(), _a.size(), _b.c_str(), _b.size(), k);
}

int lev_dist_edlib(const char *_a, int m, const char *_b, int n, int k) {
	EdlibEqualityPair *eep = NULL;
	EdlibAlignConfig config = edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, eep, 0);
	EdlibAlignResult result = edlibAlign(_a, m, _b, n, config);
	int dist = result.editDistance;
	edlibFreeAlignResult(result);
	return dist;
}