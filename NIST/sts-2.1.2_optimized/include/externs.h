
#include "../include/defs.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                   G L O B A L   D A T A  S T R U C T U R E S 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

extern BitSequence	*epsilon;				// BIT STREAM
extern BitSequence  *byteEpsilon;           // BIT STREAM IN BYTES
extern TP			tp;						// TEST PARAMETER STRUCTURE
extern FILE			*stats[NUMOFTESTS+1];	// FILE OUTPUT STREAM
extern FILE			*results[NUMOFTESTS+1];	// FILE OUTPUT STREAM
extern FILE			*freqfp;				// FILE OUTPUT STREAM
extern FILE			*summary;				// FILE OUTPUT STREAM
extern int			testVector[NUMOFTESTS+1];

extern char	generatorDir[NUMOFGENERATORS][20];
extern char	testNames[NUMOFTESTS+1][32];

//汉明距离
__attribute__((unused))
static short int LU_byte_weight[256] =
{ 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5
,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6
,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6
,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7
,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6
,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7
,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7
,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8 };

__attribute__((unused))
static unsigned int bitMask[32] = {
	0x00000000, 0x80000000, 0xC0000000, 0xE0000000,
	0xF0000000, 0xF8000000, 0xFC000000, 0xFE000000,
	0xFF000000, 0xFF800000, 0xFFC00000, 0xFFE00000,
	0xFFF00000, 0xFFF80000, 0xFFFC0000, 0xFFFE0000,
	0xFFFF0000, 0xFFFF8000, 0xFFFFC000, 0xFFFFE000,
	0xFFFFF000, 0xFFFFF800, 0xFFFFFC00, 0xFFFFFE00,
	0xFFFFFF00, 0xFFFFFF80, 0xFFFFFFC0, 0xFFFFFFE0,
	0xFFFFFFF0, 0xFFFFFFF8, 0xFFFFFFFC, 0xFFFFFFFE,
};