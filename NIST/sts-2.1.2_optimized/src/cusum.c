#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/cephes.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		    C U M U L A T I V E  S U M S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int max_plus[256] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,1,2,1,1,1,2,2,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,1,2,1,1,1,2,2,2,3,4,1,1,1,1,1,1,1,2,1,1,1,2,2,2,3,4,2,2,2,2,2,2,3,4,3,3,3,4,4,4,5,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,2,2,2,3,4,1,1,1,1,1,1,1,2,1,1,1,2,2,2,3,4,2,2,2,2,2,2,3,4,3,3,3,4,4,4,5,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,4,2,2,2,2,2,2,3,4,3,3,3,4,4,4,5,6,3,3,3,3,3,3,3,4,3,3,3,4,4,4,5,6,4,4,4,4,4,4,5,6,5,5,5,6,6,6,7,8 };
int max_minus[256] = { -8,-7,-6,-6,-6,-5,-5,-5,-6,-5,-4,-4,-4,-4,-4,-4,-6,-5,-4,-4,-4,-3,-3,-3,-4,-3,-3,-3,-3,-3,-3,-3,-6,-5,-4,-4,-4,-3,-3,-3,-4,-3,-2,-2,-2,-2,-2,-2,-4,-3,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-6,-5,-4,-4,-4,-3,-3,-3,-4,-3,-2,-2,-2,-2,-2,-2,-4,-3,-2,-2,-2,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-4,-3,-2,-2,-2,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-6,-5,-4,-4,-4,-3,-3,-3,-4,-3,-2,-2,-2,-2,-2,-2,-4,-3,-2,-2,-2,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-4,-3,-2,-2,-2,-1,-1,-1,-2,-1,0,0,0,0,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4,-3,-2,-2,-2,-1,-1,-1,-2,-1,0,0,0,0,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
int byte_sum[256] = { -8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8 };

void
CumulativeSums(int n)
{
	int		S, sup, inf, z, zrev, k;
	double	sum1, sum2, p_value;

	S = 0;
	sup = 0;
	inf = 0;
	for ( k=0; k<n; k++ ) {
		epsilon[k] ? S++ : S--;
		if ( S > sup )
			sup++;
		if ( S < inf )
			inf--;
		z = (sup > -inf) ? sup : -inf;
		zrev = (sup-S > S-inf) ? sup-S : S-inf;
	}
	
	// forward
	sum1 = 0.0;
	for ( k=(-n/z+1)/4; k<=(n/z-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*z)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*z)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/z-3)/4; k<=(n/z-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*z)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*z)/sqrt(n));
	}

	p_value = 1.0 - sum1 + sum2;
	
	fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (FORWARD) TEST\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", z);
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

	fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_CUSUM], "%f\n", p_value);
		
	// backwards
	sum1 = 0.0;
	for ( k=(-n/zrev+1)/4; k<=(n/zrev-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*zrev)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*zrev)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/zrev-3)/4; k<=(n/zrev-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*zrev)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*zrev)/sqrt(n));
	}
	p_value = 1.0 - sum1 + sum2;

	fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (REVERSE) TEST\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", zrev);
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

	fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_CUSUM]);
	fprintf(results[TEST_CUSUM], "%f\n", p_value); fflush(results[TEST_CUSUM]);
}

void
CumulativeSums1(int n)
{
	int		S, sup, inf, z, zrev, k;
	double	sum1, sum2, p_value;

	S = 0;
	sup = 0;
	inf = 0;
	int byteNum = n/8;
	for (int i = 0; i < byteNum; i++) {
		if (S + max_plus[byteEpsilon[i]] > sup) sup = S + max_plus[byteEpsilon[i]];
		if (S + max_minus[byteEpsilon[i]] < inf) inf = S + max_minus[byteEpsilon[i]];
		S += byte_sum[byteEpsilon[i]];
	}
	if (n % 8 != 0)
	{
		unsigned char local = byteEpsilon[byteNum];
		for (int j = 0; j < n - 8 * byteNum; j++) {
			((local >> (7 - j)) & 1) ? S++ : S--;
			if (S > sup)
				sup++;
			if (S < inf)
				inf--;
		}		
	}
	z = (sup > -inf) ? sup : -inf;
	zrev = (sup - S > S - inf) ? sup - S : S - inf;

	// forward
	sum1 = 0.0;
	for ( k=(-n/z+1)/4; k<=(n/z-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*z)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*z)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/z-3)/4; k<=(n/z-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*z)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*z)/sqrt(n));
	}

	p_value = 1.0 - sum1 + sum2;
	
	fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (FORWARD) TEST\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", z);
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

	fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_CUSUM], "%f\n", p_value);
		
	// backwards
	sum1 = 0.0;
	for ( k=(-n/zrev+1)/4; k<=(n/zrev-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*zrev)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*zrev)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/zrev-3)/4; k<=(n/zrev-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*zrev)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*zrev)/sqrt(n));
	}
	p_value = 1.0 - sum1 + sum2;

	fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (REVERSE) TEST\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", zrev);
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

	fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_CUSUM]);
	fprintf(results[TEST_CUSUM], "%f\n", p_value); fflush(results[TEST_CUSUM]);
}