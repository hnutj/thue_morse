#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "../include/externs.h"
#include "../include/cephes.h"

#define min(a,b) ((a) < (b) ? (a) : (b))

void
Oddhopsum(int n)
{
	int       i, ii, j, d, N, L, m, N_, K = 6, Odd;
	double    p_value,  nu[7], chi2;
    int T_;
    int fenweidian[6] = {   117,
                            121,
                            124,
                            127,
                            130,
                            134 };
    double pi[7] = { 0.1255824653,
                     0.1345129217,
                     0.1353662875,
                     0.1493542814,
                     0.1431508121,
                     0.1517626047,
                     0.1602706273 };

	BitSequence  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
	const int M = 500;
	N = (int)floor(n/M);
	if ( ((B_ = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ||
		 ((C  = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ||
		 ((P  = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ||
		 ((T  = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ) {
		printf("Insufficient Memory for Work Space:: Odd Hop Sum Test\n");
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


	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "\tODD HOP SUM\n");
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_ODDHOPSUM], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "\tNote: %d bits were discarded!\n", n%M);

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
        Odd = 1;
        T_ = 0;
		while ( N_ < M ) {
			d = (int)epsilon[ii*M+N_];
			for ( i=1; i<=L; i++ )
				d += C[i] * epsilon[ii*M+N_-i];
			d = d&1;
			if ( d == 1 ) {
				for ( i=0; i<M + 1; i++ ) {
					T[i] = C[i];
					P[i] = 0;
				}
				for ( j=0; j<M+1; j++ )
					if ( B_[j] == 1 )
						P[j+N_-m] = 1;
				for ( i=0; i<M+1; i++ )
					C[i] = (C[i] + P[i])%2;
				if ( L <= N_/2 ) {
                    if (Odd == 1) {
                        T_ += (N_ + 1 - L - L);
                    }
                    Odd = 1 - Odd;
					L = N_ + 1 - L;
					m = N_;
					for ( i=0; i<M+1; i++ )
						B_[i] = T[i];
				}
			}
			N_++;
		}
        
        i = 0;
        while (i<6 && T_ >= fenweidian[i])
            i++;
        nu[i] += 1.0;
	}
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ )
		fprintf(stats[TEST_ODDHOPSUM], "%4d ", (int)nu[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_ODDHOPSUM], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_ODDHOPSUM]);
	fprintf(results[TEST_ODDHOPSUM], "%f\n", p_value); fflush(results[TEST_ODDHOPSUM]);

	free(B_);
	free(P);
	free(C);
	free(T);
}

void
Evenhopsum(int n)
{
    int       i, ii, j, d, N, L, m, N_, K = 6, Even;
    double    p_value, nu[7], chi2;
    int T_;
    int fenweidian[6] = {   116,
                            120,
                            123,
                            126,
                            129,
                            133 };
    double pi[7] = {    0.1255824653,
                        0.1345129217,
                        0.1353662875,
                        0.1493542814,
                        0.1431508121,
                        0.1517626047,
                        0.1602706273 };

    BitSequence  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
    const int M = 500;
    N = (int)floor(n / M);
    if (((B_ = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
        ((C = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
        ((P = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
        ((T = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL)) {
        printf("Insufficient Memory for Work Space:: Even Hop Sum Test\n");
        if (B_ != NULL)
            free(B_);
        if (C != NULL)
            free(C);
        if (P != NULL)
            free(P);
        if (T != NULL)
            free(T);
        return;
    }


    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "\tEVEN HOP SUM\n");
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_EVENHOPSUM], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "        F R E Q U E N C Y                            \n");
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "\tNote: %d bits were discarded!\n", n%M);

    for (i = 0; i<K + 1; i++)
        nu[i] = 0.00;
    for (ii = 0; ii<N; ii++) {
        for (i = 0; i<M; i++) {
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
        Even = 0;
        T_ = 0;
        while (N_ < M) {
            d = (int)epsilon[ii*M + N_];
            for (i = 1; i <= L; i++)
                d += C[i] * epsilon[ii*M + N_ - i];
            d = d & 1;
            if (d == 1) {
                for (i = 0; i<M + 1; i++) {
                    T[i] = C[i];
                    P[i] = 0;
                }
                for (j = 0; j<M + 1; j++)
                    if (B_[j] == 1)
                        P[j + N_ - m] = 1;
                for (i = 0; i<M + 1; i++)
                    C[i] = (C[i] + P[i]) % 2;
                if (L <= N_ / 2) {
                    if (Even == 1) {
                        T_ += (N_ + 1 - L - L);
                    }
                    Even = 1 - Even;
                    L = N_ + 1 - L;
                    m = N_;
                    for (i = 0; i<M + 1; i++)
                        B_[i] = T[i];
                }
            }
            N_++;
        }

        i = 0;
        while (i<6 && T_ >= fenweidian[i])
            i++;
        nu[i] += 1.0;
    }
    chi2 = 0.00;
    for (i = 0; i<K + 1; i++)
        fprintf(stats[TEST_EVENHOPSUM], "%4d ", (int)nu[i]);
    for (i = 0; i<K + 1; i++)
        chi2 += pow(nu[i] - N*pi[i], 2) / (N*pi[i]);
    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);

    fprintf(stats[TEST_EVENHOPSUM], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_EVENHOPSUM]);
    fprintf(results[TEST_EVENHOPSUM], "%f\n", p_value); fflush(results[TEST_EVENHOPSUM]);

    free(B_);
    free(P);
    free(C);
    free(T);
}

void
Jump(int n)
{
    int       i, ii, j, d, N, L, m, N_, K = 6;
    double    p_value, nu[7], chi2;
    int T_;
    int fenweidian[6] = {   117,
                            121,
                            124,
                            127,
                            130,
                            134 };
    double pi[7] = {    0.1314787313,
                        0.1386850066,
                        0.1380215051,
                        0.1505690987,
                        0.1424525791,
                        0.1484288066,
                        0.1503642727 };

    BitSequence  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
    const int M = 500;
    N = (int)floor(n / M);
    if (((B_ = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
    ((C = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
    ((P = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
    ((T = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL)) {
        printf("Insufficient Memory for Work Space:: Jump Complexity Test\n");
        if (B_ != NULL)
            free(B_);
        if (C != NULL)
            free(C);
        if (P != NULL)
            free(P);
        if (T != NULL)
            free(T);
        return;
    }


    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "\tJUMP\n");
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_JUMP], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "        F R E Q U E N C Y                            \n");
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "\tNote: %d bits were discarded!\n", n%M);

    for (i = 0; i<K + 1; i++)
        nu[i] = 0.00;
    for (ii = 0; ii<N; ii++) {
        for (i = 0; i<M; i++) {
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
        T_ = 0;
        while (N_ < M) {
            d = (int)epsilon[ii*M + N_];
            for (i = 1; i <= L; i++)
                d += C[i] * epsilon[ii*M + N_ - i];
            d = d & 1;
            if (d == 1) {
                for (i = 0; i<M + 1; i++) {
                    T[i] = C[i];
                    P[i] = 0;
                }
                for (j = 0; j<M + 1; j++)
                    if (B_[j] == 1)
                        P[j + N_ - m] = 1;
                for (i = 0; i<M + 1; i++)
                    C[i] = (C[i] + P[i]) % 2;
                
                // 当线性复杂度L变化的时候，就有一个jump
                if (L <= N_ / 2) {
                    
                    T_++;
                    L = N_ + 1 - L;
                    m = N_;
                    for (i = 0; i<M + 1; i++)
                        B_[i] = T[i];
                }
            }
            N_++;
        }

        i = 0;
        while (i<6 && T_ >= fenweidian[i])
            i++;
        nu[i] += 1.0;
    }
    chi2 = 0.00;
    for (i = 0; i<K + 1; i++)
        fprintf(stats[TEST_JUMP], "%4d ", (int)nu[i]);
    for (i = 0; i<K + 1; i++)
        chi2 += pow(nu[i] - N*pi[i], 2) / (N*pi[i]);
    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);

    fprintf(stats[TEST_JUMP], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_JUMP]);
    fprintf(results[TEST_JUMP], "%f\n", p_value); fflush(results[TEST_JUMP]);

    free(B_);
    free(P);
    free(C);
    free(T);
}

void
Oddhopsum1(int n)
{
	int       i, ii, N, K = 6;
	double    p_value,  nu[7], chi2;
    int fenweidian[6] = {   117,
                            121,
                            124,
                            127,
                            130,
                            134 };
    double pi[7] = { 0.1255824653,
                     0.1345129217,
                     0.1353662875,
                     0.1493542814,
                     0.1431508121,
                     0.1517626047,
                     0.1602706273 };

	int 	  *counter;
	unsigned int    power2[32], lenS, bitLenS;
	const int M = 500;
	N = (int)floor(n/M);
	lenS = M; 
	bitLenS = (lenS + 31) / 32;
	power2[0] = 1;
	for (int i = 1; i < 32; i++)
		power2[i] = 2 * power2[i - 1];
	if( (counter = (int*)malloc(sizeof(int) * N)) == NULL ){
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		return;		
	}

	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "\tODD HOP SUM\n");
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_ODDHOPSUM], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "\tNote: %d bits were discarded!\n", n%M);

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
		int lenC, lenD, m, n, q, r, d, upperBound, wordCnt, shiftD, startC, T_, Odd;
		// int word, bitPos;
		n = lenC = lenD = 0;
		m = -1;
        Odd = 1;
        T_ = 0;
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
                    if (Odd == 1) {
                        // T_+=L_new-L_old=(n+1-lenC)-lenC
                        T_ += n + 1 - lenC - lenC;
                    }
                    Odd = 1 - Odd;
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
		counter[ii] = T_;

		free(bitC);
		free(--bitD);
		free(bitTmp);

		for (int i = 0; i < 32; i++)
			free(bitRS[i]);
	}

	for ( i=0; i<K+1; i++ )
		nu[i] = 0.00;
    for( ii=0; ii<N; ii++){
        i = 0;
        int T_ = counter[ii];
        while (i<6 && T_ >= fenweidian[i])
            i++;
        nu[i] += 1.0;
	}
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ )
		fprintf(stats[TEST_ODDHOPSUM], "%4d ", (int)nu[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_ODDHOPSUM], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_ODDHOPSUM]);
	fprintf(results[TEST_ODDHOPSUM], "%f\n", p_value); fflush(results[TEST_ODDHOPSUM]);

}

void
Evenhopsum1(int n)
{
	int       i, ii, N, K = 6;
    double    p_value, nu[7], chi2;
    int fenweidian[6] = {   116,
                            120,
                            123,
                            126,
                            129,
                            133 };
    double pi[7] = {    0.1255824653,
                        0.1345129217,
                        0.1353662875,
                        0.1493542814,
                        0.1431508121,
                        0.1517626047,
                        0.1602706273 };
  
    const int M = 500;
    N = (int)floor(n / M);
  	int 	  *counter;
	unsigned int    power2[32], lenS, bitLenS;
	lenS = M; 
	bitLenS = (lenS + 31) / 32;
	power2[0] = 1;
	for (int i = 1; i < 32; i++)
		power2[i] = 2 * power2[i - 1];
	if( (counter = (int*)malloc(sizeof(int) * N)) == NULL ){
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		return;		
	}

    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "\tEVEN HOP SUM\n");
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_EVENHOPSUM], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "        F R E Q U E N C Y                            \n");
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "\tNote: %d bits were discarded!\n", n%M);

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
		int lenC, lenD, m, n, q, r, d, upperBound, wordCnt, shiftD, startC, T_, Even;
		// int word, bitPos;
		n = lenC = lenD = 0;
		m = -1;
        Even = 0;
        T_ = 0;
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
                    if (Even == 1) {
                        T_ += (n + 1 - lenC - lenC);
                    }
                    Even = 1 - Even;
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
		counter[ii] = T_;

		free(bitC);
		free(--bitD);
		free(bitTmp);

		for (int i = 0; i < 32; i++)
			free(bitRS[i]);
	}

    for (i = 0; i<K + 1; i++)
        nu[i] = 0.00;
    for (ii = 0; ii<N; ii++) {
        int T_=counter[ii];
        i = 0;
        while (i<6 && T_ >= fenweidian[i])
            i++;
        nu[i] += 1.0;
    }
    chi2 = 0.00;
    for (i = 0; i<K + 1; i++)
        fprintf(stats[TEST_EVENHOPSUM], "%4d ", (int)nu[i]);
    for (i = 0; i<K + 1; i++)
        chi2 += pow(nu[i] - N*pi[i], 2) / (N*pi[i]);
    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);

    fprintf(stats[TEST_EVENHOPSUM], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_EVENHOPSUM]);
    fprintf(results[TEST_EVENHOPSUM], "%f\n", p_value); fflush(results[TEST_EVENHOPSUM]);

}

void
Jump1(int n)
{
	int       i, ii, N, K = 6;
    double    p_value, nu[7], chi2;
    int 	  *counter;
	unsigned int    power2[32], lenS, bitLenS;

    int fenweidian[6] = {   117,
                            121,
                            124,
                            127,
                            130,
                            134 };
    double pi[7] = {    0.1314787313,
                        0.1386850066,
                        0.1380215051,
                        0.1505690987,
                        0.1424525791,
                        0.1484288066,
                        0.1503642727 };

    const int M = 500;
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

    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "\tJUMP\n");
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_JUMP], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "        F R E Q U E N C Y                            \n");
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "\tNote: %d bits were discarded!\n", n%M);

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
		int lenC, lenD, m, n, q, r, d, upperBound, wordCnt, shiftD, startC, T_;
		// int word, bitPos;
		n = lenC = lenD = 0;
		m = -1;
        T_ = 0;
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

				//更新D和m,每当线性复杂度lenC变化，跳跃T_++
				if (lenC <= (n >> 1)) {
                    T_++;
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
		counter[ii] = T_;

		free(bitC);
		free(--bitD);
		free(bitTmp);

		for (int i = 0; i < 32; i++)
			free(bitRS[i]);
	}

    for (i = 0; i<K + 1; i++)
        nu[i] = 0.00;
    for (ii = 0; ii<N; ii++) {
        int T_ = counter[ii];
        i = 0;
        while (i<6 && T_ >= fenweidian[i])
            i++;
        nu[i] += 1.0;
    }

    chi2 = 0.00;
    for (i = 0; i<K + 1; i++)
        fprintf(stats[TEST_JUMP], "%4d ", (int)nu[i]);
    for (i = 0; i<K + 1; i++)
        chi2 += pow(nu[i] - N*pi[i], 2) / (N*pi[i]);
    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);

    fprintf(stats[TEST_JUMP], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_JUMP]);
    fprintf(results[TEST_JUMP], "%f\n", p_value); fflush(results[TEST_JUMP]);

}

// jump / OddHop / EvenHop
void
JumpHopComplexity1(int n)
{
	int       i, ii, N, K = 6;
	double    p_value,  nu[7], chi2;
	int 	  *counter[3];
	unsigned int    power2[32], lenS, bitLenS;
	const int       M = 500;

    int oddHop_fenweidian[6] = { 117, 121, 124, 127, 130, 134 }; //in line with jump test
    int evenHop_fenweidian[6] = { 116, 120, 123, 126, 129, 133 };
    double hop_pi[7] = { 0.1255824653, 0.1345129217, 0.1353662875,
                     0.1493542814, 0.1431508121, 0.1517626047, 0.1602706273 };
    double jump_pi[7] = { 0.1314787313, 0.1386850066, 0.1380215051,
                        0.1505690987, 0.1424525791, 0.1484288066, 0.1503642727 };
	N = (int)floor(n/M);
	lenS = M; 
	bitLenS = (lenS + 31) / 32;
	power2[0] = 1;
	for (int i = 1; i < 32; i++)
		power2[i] = 2 * power2[i - 1];
    for( i = 0; i < 3; i++){
        if( (counter[i] = (int*)malloc(sizeof(int) * N)) == NULL ){
            printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
            return;		
        }        
    }

	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "\tODD HOP SUM\n");
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_ODDHOPSUM], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
	fprintf(stats[TEST_ODDHOPSUM], "-----------------------------------------------------\n");
	fprintf(stats[TEST_ODDHOPSUM], "\tNote: %d bits were discarded!\n", n%M);

    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "\tEVEN HOP SUM\n");
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_EVENHOPSUM], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "        F R E Q U E N C Y                            \n");
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
    fprintf(stats[TEST_EVENHOPSUM], "-----------------------------------------------------\n");
    fprintf(stats[TEST_EVENHOPSUM], "\tNote: %d bits were discarded!\n", n%M);

    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "\tJUMP\n");
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_JUMP], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "        F R E Q U E N C Y                            \n");
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
    fprintf(stats[TEST_JUMP], "-----------------------------------------------------\n");
    fprintf(stats[TEST_JUMP], "\tNote: %d bits were discarded!\n", n%M);

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
		int lenC, lenD, m, n, q, r, d, upperBound, wordCnt, shiftD, startC, T_odd, T_even, T_jump, Odd;
		// int word, bitPos;
		n = lenC = lenD = 0;
		m = -1;
        Odd = 1;
        T_odd = 0;
        T_even = 0;
        T_jump = 0;
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
                    T_jump++;
                    if (Odd == 1) {
                        // T_odd+=L_new-L_old=(n+1-lenC)-lenC
                        T_odd += n + 1 - lenC - lenC;
                    }else{
                        T_even += n + 1 - lenC - lenC;
                    }
                    Odd = 1 - Odd;
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
		counter[0][ii] = T_odd;
		counter[1][ii] = T_even;
        counter[2][ii] = T_jump;

		free(bitC);
		free(--bitD);
		free(bitTmp);

		for (int i = 0; i < 32; i++)
			free(bitRS[i]);
	}

    // oddHop test
	for ( i=0; i<K+1; i++ )
		nu[i] = 0.00;
    for( ii=0; ii<N; ii++){
        i = 0;
        int T_ = counter[0][ii];
        while (i<6 && T_ >= oddHop_fenweidian[i])
            i++;
        nu[i] += 1.0;
	}
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ )
		fprintf(stats[TEST_ODDHOPSUM], "%4d ", (int)nu[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*hop_pi[i], 2) / (N*hop_pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_ODDHOPSUM], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_ODDHOPSUM]);
	fprintf(results[TEST_ODDHOPSUM], "%f\n", p_value); fflush(results[TEST_ODDHOPSUM]);

    // evenHop test
    for (i = 0; i<K + 1; i++)
        nu[i] = 0.00;
    for (ii = 0; ii<N; ii++) {
        int T_=counter[1][ii];
        i = 0;
        while (i<6 && T_ >= evenHop_fenweidian[i])
            i++;
        nu[i] += 1.0;
    }
    chi2 = 0.00;
    for (i = 0; i<K + 1; i++)
        fprintf(stats[TEST_EVENHOPSUM], "%4d ", (int)nu[i]);
    for (i = 0; i<K + 1; i++)
        chi2 += pow(nu[i] - N*hop_pi[i], 2) / (N*hop_pi[i]);
    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);

    fprintf(stats[TEST_EVENHOPSUM], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_EVENHOPSUM]);
    fprintf(results[TEST_EVENHOPSUM], "%f\n", p_value); fflush(results[TEST_EVENHOPSUM]);

    // jump test
    for (i = 0; i<K + 1; i++)
        nu[i] = 0.00;
    for (ii = 0; ii<N; ii++) {
        int T_ = counter[2][ii];
        i = 0;
        while (i<6 && T_ >= oddHop_fenweidian[i])
            i++;
        nu[i] += 1.0;
    }
    chi2 = 0.00;
    for (i = 0; i<K + 1; i++)
        fprintf(stats[TEST_JUMP], "%4d ", (int)nu[i]);
    for (i = 0; i<K + 1; i++)
        chi2 += pow(nu[i] - N*jump_pi[i], 2) / (N*jump_pi[i]);
    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);

    fprintf(stats[TEST_JUMP], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_JUMP]);
    fprintf(results[TEST_JUMP], "%f\n", p_value); fflush(results[TEST_JUMP]);

}