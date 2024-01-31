#ifndef __constants_h
#define __constants_h

#include <climits>

// Represents an unknown (or unimportant) distance between blocks.
#define INF_DIST INT_MAX

// The cost of a character deletion operation. Currently set to 1.
#define C_CHAR_DELETION 1
// The cost of a character insertion operation. Currently set to 1.
#define C_CHAR_INSERTION 1
// The cost of a single nucleotide polymorphism (SNP) operation. Currently set to 1.
#define C_CHAR_SNP 1

// The cost of a block deletion operation. Currently set to 1.
#define C_BLOCK_DELETION 1
// The cost of a block move operation. Currently set to 1.
#define C_BLOCK_MOVE 1
// The cost of a block reverse operation. Currently set to 1.
#define C_BLOCK_REVERSE 1

// The default minimum block size for the alignment. Currently set to 8.
#define DEFAULT_L_MIN 8
// The default maximum block size for the alignment. Currently set to 15.
#define DEFAULT_L_MAX 15
// The default maximum number of iterations for the alignment algorithm. Currently set to 3.
#define DEFAULT_MAX_ITERATIONS 3
// The default maximum error rate for the algorithm. Currently set to 0.3.
#define DEFAULT_ERROR_RATE 0.3



#endif //__constants_h