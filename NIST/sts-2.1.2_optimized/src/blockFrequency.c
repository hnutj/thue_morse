#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/cephes.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                    B L O C K  F R E Q U E N C Y  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
BlockFrequency(int M, int n)
{
	int		i, j, N, blockSum;
	double	p_value, sum, pi, v, chi_squared;
	
	N = n/M; 		/* # OF SUBSTRING BLOCKS      */
	sum = 0.0;
	
	for ( i=0; i<N; i++ ) {
		blockSum = 0;
		for ( j=0; j<M; j++ )
			blockSum += epsilon[j+i*M];
		pi = (double)blockSum/(double)M;
		v = pi - 0.5;
		sum += v*v;
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N/2.0, chi_squared/2.0);

	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t\tBLOCK FREQUENCY TEST\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(a) Chi^2           = %f\n", chi_squared);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(b) # of substrings = %d\n", N);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(c) block length    = %d\n", M);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(d) Note: %d bits were discarded.\n", n % M);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_BLOCK_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BLOCK_FREQUENCY]);
	fprintf(results[TEST_BLOCK_FREQUENCY], "%f\n", p_value); fflush(results[TEST_BLOCK_FREQUENCY]);
}

// 块内频率检测
void
BlockFrequency1(int M, int n)
{
	int		i, j, k, N, blockSum;
	double	p_value, sum, pi, v, chi_squared;
	
	N = n/M; 		/* # OF SUBSTRING BLOCKS      */
	sum = 0.0;
	
	for ( i=0; i<N; i++ ) {
		blockSum = 0;
		// for ( j=0; j<M; j++ )
		// 	blockSum += epsilon[j+i*M];
		int beginPos = i * M, endPos = (i + 1) * M - 1;
		int byteBeginPos = beginPos / 8, byteEndPos = endPos / 8;
		for (j = byteBeginPos; j < byteEndPos; j++) {
			blockSum += LU_byte_weight[byteEpsilon[j]];
		}
		//减去前面多余的
		int left1 = beginPos - 8 * byteBeginPos;
		for (k = 0; k < left1; k++) {
			blockSum -= ((byteEpsilon[byteBeginPos] >> (7 - k)) & 1);
		}
		//加上后面不足的
		int left2 = endPos - byteEndPos * 8 + 1;
		for (k = 0; k < left2; k++) {
			blockSum += ((byteEpsilon[byteEndPos] >> (7 - k)) & 1);
		}
		pi = (double)blockSum/(double)M;
		v = pi - 0.5;
		sum += v*v;
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N/2.0, chi_squared/2.0);

	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t\tBLOCK FREQUENCY TEST\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(a) Chi^2           = %f\n", chi_squared);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(b) # of substrings = %d\n", N);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(c) block length    = %d\n", M);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(d) Note: %d bits were discarded.\n", n % M);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_BLOCK_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BLOCK_FREQUENCY]);
	fprintf(results[TEST_BLOCK_FREQUENCY], "%f\n", p_value); fflush(results[TEST_BLOCK_FREQUENCY]);
}
