#ifndef __constants_h
#define __constants_h

#include <climits>

// Unknown distance between between blocks
#define INF_DIST INT_MAX
#define UNKNOWN_TUPLE std::make_tuple(-1,-1, false)

// Character Operation Costs
#define C_CHAR_DELETION 1
#define C_CHAR_INSERTION 1
#define C_CHAR_SNP 1

// Block Operation Costs
#define C_BLOCK_DELETION 1
#define C_BLOCK_MOVE 1
#define C_BLOCK_REVERSE 1

#endif //__constants_h