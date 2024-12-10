#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "../include/externs.h"
#include "../include/cephes.h"  

#define min(a,b) ((a) < (b) ? (a) : (b))

void
LinearComplexity(int M, int n)
{
	int       i, ii, j, d, N, L, m, N_, parity, sign, K = 6;
	double    p_value, T_, mean, nu[7], chi2;
	double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
	BitSequence  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
	
	N = (int)floor(n/M);
	if ( ((B_ = (BitSequence *) calloc(M, sizeof(BitSequence))) == NULL) ||
		 ((C  = (BitSequence *) calloc(M, sizeof(BitSequence))) == NULL) ||
		 ((P  = (BitSequence *) calloc(M, sizeof(BitSequence))) == NULL) ||
		 ((T  = (BitSequence *) calloc(M, sizeof(BitSequence))) == NULL) ) {
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		if ( B_ != NULL )
			free(B_);
		if ( C != NULL )
			free(C);
		if ( P != NULL )
			free(P);
		if ( T != NULL )
			free(T);
		return;
	}


	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tL I N E A R  C O M P L E X I T Y\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "  C0   C1   C2   C3   C4   C5   C6    CHI2    P-value\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tNote: %d bits were discarded!\n", n%M);

	for ( i=0; i<K+1; i++ )
		nu[i] = 0.00;
	for ( ii=0; ii<N; ii++ ) {
		for ( i=0; i<M; i++ ) {
			B_[i] = 0;
			C[i] = 0;
			T[i] = 0;
			P[i] = 0;
		}
		L = 0;
		m = -1;
		d = 0;
		C[0] = 1;
		B_[0] = 1;
		
		/* DETERMINE LINEAR COMPLEXITY */
		N_ = 0;
		while ( N_ < M ) {
			d = (int)epsilon[ii*M+N_];
			for ( i=1; i<=L; i++ )
				d += C[i] * epsilon[ii*M+N_-i];
			d = d%2;
			if ( d == 1 ) {
				for ( i=0; i<M; i++ ) {
					T[i] = C[i];
					P[i] = 0;
				}
				for ( j=0; j<M; j++ )
					if ( B_[j] == 1 )
						P[j+N_-m] = 1;
				for ( i=0; i<M; i++ )
					C[i] = (C[i] + P[i])%2;
				if ( L <= N_/2 ) {
					L = N_ + 1 - L;
					m = N_;
					for ( i=0; i<M; i++ )
						B_[i] = T[i];
				}
			}
			N_++;
		}
		if ( (parity = (M+1)%2) == 0 ) 
			sign = -1;
		else 
			sign = 1;
		mean = M/2.0 + (9.0+sign)/36.0 - 1.0/pow(2, M) * (M/3.0 + 2.0/9.0);
		if ( (parity = M%2) == 0 )
			sign = 1;
		else 
			sign = -1;
		T_ = sign * (L - mean) + 2.0/9.0;
		
		if ( T_ <= -2.5 )
			nu[0]++;
		else if ( T_ > -2.5 && T_ <= -1.5 )
			nu[1]++;
		else if ( T_ > -1.5 && T_ <= -0.5 )
			nu[2]++;
		else if ( T_ > -0.5 && T_ <= 0.5 )
			nu[3]++;
		else if ( T_ > 0.5 && T_ <= 1.5 )
			nu[4]++;
		else if ( T_ > 1.5 && T_ <= 2.5 )
			nu[5]++;
		else
			nu[6]++;
	}
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ ) 
		fprintf(stats[TEST_LINEARCOMPLEXITY], "%4d ", (int)nu[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_LINEARCOMPLEXITY], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_LINEARCOMPLEXITY]);
	fprintf(results[TEST_LINEARCOMPLEXITY], "%f\n", p_value); fflush(results[TEST_LINEARCOMPLEXITY]);

	free(B_);
	free(P);
	free(C);
	free(T);
}

void
LinearComplexity1(int M, int n)
{
	int       i, N, L, parity, sign, K = 6;
	double    p_value, T_, mean, nu[7], chi2;
	//校正
	//double    pi[7] = { 0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
	double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
	int 	  *counter;
	unsigned int    power2[32], lenS, bitLenS;

	N = (int)floor(n/M); //分块块数
	lenS = M; 
	bitLenS = (lenS + 31) / 32;
	power2[0] = 1;
	for (int i = 1; i < 32; i++)
		power2[i] = 2 * power2[i - 1];
	if( (counter = (int*)malloc(sizeof(int) * N)) == NULL ){
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		return;		
	}

	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tL I N E A R  C O M P L E X I T Y\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "  C0   C1   C2   C3   C4   C5   C6    CHI2    P-value\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tNote: %d bits were discarded!\n", n%M);
	
	#pragma omp parallel for 
	for (int ii = 0; ii < N; ii++) {

		unsigned int* bitRS[32], * bitC, * bitD, * bitTmp;
		for (int i = 0; i < 32; i++)
			bitRS[i] = (unsigned int*)malloc(sizeof(unsigned int) * (bitLenS + 1));
		bitTmp = (unsigned int*)malloc(sizeof(unsigned int) * bitLenS); 
		bitC = (unsigned int*)malloc(sizeof(unsigned int) * bitLenS);
		bitD = (unsigned int*)malloc(sizeof(unsigned int) * (bitLenS + 1));
		bitD[0] = 0;
		bitD++;

		if ( (bitRS == NULL) || (bitTmp == NULL) || (bitC == NULL) || (bitD == NULL )) {
			printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		}

		// 将这个块中对应的序列存入bitRS，初始化bitC和bitD
		int endPos = (ii + 1) * M - 1;
		bitD[0] = bitC[0] = power2[31];
		for (int i = 1; i < bitLenS; i++)
			bitD[i] = bitC[i] = 0;

		int j = 0;
		bitRS[0][j] = 0;
		int bitPos = 31;
		// bitRS[0]以逆序方式存储每个子块的序列，越靠近endPos存储在int的高位
		for (int i = 0; i < lenS; i++) {
			if (epsilon[endPos - i])
				bitRS[0][j] |= power2[bitPos];
			bitPos--;
			if (bitPos == -1) {
				bitPos = 31;
				j++;
				bitRS[0][j] = 0;
			}
		}
		bitRS[0][bitLenS] = 0;
		// bitRS[i]以逆序方式存储beginPos~endPos-i这一串序列的值
		for (int i = 1; i < 32; i++) {
			for (j = 0; j < bitLenS; j++)
				bitRS[i][j] = (bitRS[0][j] << i) |
				((bitRS[0][j + 1] & bitMask[i]) >> (32 - i));
			bitRS[i][bitLenS] = 0;
		}

		// C存储最终的特征多项式，l=lenC表示特征多项式的degree/多项式C的次数
		// D存储中间临时的特征多项式，lenD为次数，m代表S序列的前m个数可以用多项式D生成，D用以更新f_n
		int lenC, lenD, m, n, q, r, d, upperBound, wordCnt, shiftD, startC;
		// int word, bitPos;
		n = lenC = lenD = 0;
		m = -1;
		while (n < lenS) {
			q = (lenS - 1 - n) >> 5;
			r = (lenS - 1 - n) & ~bitMask[27];
			d = 0;
			for (int i = 0; i < (lenC + 32) >> 5; i++)
				d ^= bitC[i] & bitRS[r][q + i];
			d = d - ((d >> 1) & 0x55555555);
			d = (d & 0x33333333) + ((d >> 2) & 0x33333333);
			d = (((d + (d >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
			// 如果d=0，代表前n位能够用之前的特征多项式f_(n-1)表示，继续推进即可
			// 如果d=1，代表不能表示，要生成f_n这个新的特征多项式,l_(n-1)与l_n必定不同
			// 此时区分2l>n和2l<=n这两种情况，只有在后者情况，D和l=lenC才会有变化
			if (d & 1) {
				if (lenC <= (n >> 1))
					for (int i = 0; i < (lenC + 32) >> 5; i++)
						bitTmp[i] = bitC[i];

				//更新bitC
				upperBound = min(lenD + 1, (int)lenS + m - n);
				startC = (n - m) >> 5;
				shiftD = (n - m) & ~bitMask[27];
				wordCnt = 0;
				if (shiftD) {
					upperBound -= (32 - shiftD);
					wordCnt++;
				}
				wordCnt += (upperBound + 31) >> 5;
				for (int i = 0; i < wordCnt; i++)
					if (shiftD)
						bitC[startC + i] ^= ((bitD[i - 1] & ~bitMask[32 - shiftD]) << (32 - shiftD))
						| ((bitD[i] & bitMask[32 - shiftD]) >> shiftD);
					else
						bitC[startC + i] ^= bitD[i];

				//更新D和m
				if (lenC <= (n >> 1)) {
					for (int i = 0; i < (lenC + 32) >> 5; i++)
						bitD[i] = bitTmp[i];
					lenD = lenC;
					lenC = n + 1 - lenC;
					m = n;
				}
			}
			n++;
		}

		//存入
		counter[ii] = lenC;

		free(bitC);
		free(--bitD);
		free(bitTmp);

		for (int i = 0; i < 32; i++)
			free(bitRS[i]);
	}

	for ( i=0; i<K+1; i++ )
		nu[i] = 0.00;
	if ( (parity = (M+1)%2) == 0 ) 
		sign = -1;
	else 
		sign = 1;
	mean = M/2.0 + (9.0+sign)/36.0 - 1.0/pow(2, M) * (M/3.0 + 2.0/9.0);
	if ( (parity = M%2) == 0 )
		sign = 1;
	else 
		sign = -1;
	for( i=0; i<N; i++){
		L = counter[i];
		T_ = sign * (L - mean) + 2.0/9.0;
		
		if ( T_ <= -2.5 )
			nu[0]++;
		else if ( T_ > -2.5 && T_ <= -1.5 )
			nu[1]++;
		else if ( T_ > -1.5 && T_ <= -0.5 )
			nu[2]++;
		else if ( T_ > -0.5 && T_ <= 0.5 )
			nu[3]++;
		else if ( T_ > 0.5 && T_ <= 1.5 )
			nu[4]++;
		else if ( T_ > 1.5 && T_ <= 2.5 )
			nu[5]++;
		else
			nu[6]++;	
	}

	chi2 = 0.00;
	for ( i=0; i<K+1; i++ ) 
		fprintf(stats[TEST_LINEARCOMPLEXITY], "%4d ", (int)nu[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_LINEARCOMPLEXITY], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_LINEARCOMPLEXITY]);
	fprintf(results[TEST_LINEARCOMPLEXITY], "%f\n", p_value); fflush(results[TEST_LINEARCOMPLEXITY]);
}