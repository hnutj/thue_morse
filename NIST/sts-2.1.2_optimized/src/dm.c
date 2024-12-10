#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "../include/externs.h"
#include "../include/cephes.h"
/*
    DM1 using chi square to calculate the P-value
*/

#define min(a,b) ((a) < (b) ? (a) : (b))

void
DeviationMeasure1CS(int n)
{
	int       i, ii, j, d, N, L, m, N_, K = 6;
	double    p_value, T_ , nu[7], chi2;
    double fenweidian[6] = { 335.0,
                             351.0,
                             364.0,
                             377.0,
                             393.0,
                             415.0 };
    double    pi[7]={ 0.1400018365,
                      0.1422095933,
                      0.1394350127,
                      0.1403124093,
                      0.1503319371,
                      0.1444503240,
                      0.1432588871 };

	BitSequence  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
	const int M = 500;
	N = (int)floor(n/M);
	if ( ((B_ = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ||
		 ((C  = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ||
		 ((P  = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ||
		 ((T  = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ) {
		printf("Insufficient Memory for Work Space:: Deviation Measure I CS Test\n");
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


	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tDEVIATION MEASURE I CS\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tNote: %d bits were discarded!\n", n%M);

	for ( i=0; i<K+1; i++ )
		nu[i] = 0.00;
	for ( ii=0; ii<N; ii++ ) {
        T_ = 0.0;
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
					L = N_ + 1 - L;
					m = N_;
					for ( i=0; i<M+1; i++ )
						B_[i] = T[i];
				}
			}
			T_ += fabs(L - (N_+1)/2.0);
			N_++;
		}

		/*if ( T_ < fenweidian[0])
			nu[0]+=1.0;
		else if ( T_ < fenweidian[1] )
			nu[1]+=1.0;
		else if ( T_ < fenweidian[2])
			nu[2]+=1.0;
		else if ( T_ < fenweidian[3])
			nu[3]+=1.0;
        else if (T_ < fenweidian[4])
            nu[4] += 1.0;
        else if (T_ < fenweidian[5])
            nu[5] += 1.0;
		else
			nu[6]+=1.0;*/
        
        i = 0;
        while (i<6 && T_ >= fenweidian[i])
            i++;
        nu[i] += 1.0;
	}
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ )
		fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "%4d ", (int)nu[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_DEVIATIONMEARSURE1CS]);
	fprintf(results[TEST_DEVIATIONMEARSURE1CS], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE1CS]);

	free(B_);
	free(P);
	free(C);
	free(T);
}
/*
DM1 using binomial distribution to calculate the P-value
the percentile is chose to be upper 50% 
*/
void
DeviationMeasure1(int n)
{
    int       i, ii, j, d, N, L, m, N_;
    double    p_value, T_, nu;
    double    fenweidian = 371.0;
    double    pi = 0.4909290857;

    BitSequence  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
    const int M = 500;
    N = (int)floor(n / M);
    if (((B_ = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
        ((C = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
        ((P = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
        ((T = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL)) {
        printf("Insufficient Memory for Work Space:: Deviation Measure I Test\n");
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


    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tDEVIATION MEASURE I\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tNote: %d bits were discarded!\n", n%M);
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "  Frequency     CHI2      P-value\n");
    
    nu = 0.00;
    for (ii = 0; ii<N; ii++) {
        T_ = 0.0;
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
                    L = N_ + 1 - L;
                    m = N_;
                    for (i = 0; i<M + 1; i++)
                        B_[i] = T[i];
                }
            }
            T_ += fabs(L - (N_ + 1) / 2.0);
            N_++;
        }

       

        if(T_ > fenweidian)
            nu += 1.0;
    }
    double p = nu / N;
    double z = (p - pi) / sqrt(pi * (1 - pi) / N);
    
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "    %4d    ", (int)nu);
    p_value = cephes_erfc(fabs(z) / sqrt(2));

    fprintf(stats[TEST_DEVIATIONMEARSURE1], "%9.6f    %9.6f\n", z, p_value); fflush(stats[TEST_DEVIATIONMEARSURE1]);
    fprintf(results[TEST_DEVIATIONMEARSURE1], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE1]);

    free(B_);
    free(P);
    free(C);
    free(T);
}


/*
DM2 using chi square to calculate the P-value
*/
void
DeviationMeasure2CS(int n)
{
	int       i, ii, j, d, N, L, m, N_, K = 6;
	double    p_value, T_ , nu[7], chi2;
    double fenweidian[6] = {
        422.5,
        469.5,
        511.5,
        556.5,
        612.5,
        699.5
    };
	double    pi[7]={   0.1404715615,
                        0.1427569687,
                        0.1430949506,
                        0.1436457382,
                        0.1440210422,
                        0.1423812606,
                        0.1436284782 };
	BitSequence  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
	const int M = 500;
	N = (int)floor(n/M);
	if ( ((B_ = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ||
		 ((C  = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ||
		 ((P  = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ||
		 ((T  = (BitSequence *) calloc(M+1, sizeof(BitSequence))) == NULL) ) {
		printf("Insufficient Memory for Work Space:: Deviation Measure II CS Test\n");
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


	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tDEVIATION MEASURE II CS\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tNote: %d bits were discarded!\n", n%M);

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
        T_ = 0.0;
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
					L = N_ + 1 - L;
					m = N_;
					for ( i=0; i<M+1; i++ )
						B_[i] = T[i];
				}
			}
			T_ += (L - (N_+1)/2.0)*(L - (N_+1)/2.0);
			N_++;
		}

        /*if ( T_ < fenweidian[0])
        nu[0]+=1.0;
        else if ( T_ < fenweidian[1] )
        nu[1]+=1.0;
        else if ( T_ < fenweidian[2])
        nu[2]+=1.0;
        else if ( T_ < fenweidian[3])
        nu[3]+=1.0;
        else if (T_ < fenweidian[4])
        nu[4] += 1.0;
        else if (T_ < fenweidian[5])
        nu[5] += 1.0;
        else
        nu[6]+=1.0;*/

        i = 0;
        while (i<6 && T_ >= fenweidian[i])
            i++;
        nu[i] += 1.0;
	}
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ )
		fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "%4d ", (int)nu[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_DEVIATIONMEARSURE2CS]);
	fprintf(results[TEST_DEVIATIONMEARSURE2CS], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE2CS]);

	free(B_);
	free(P);
	free(C);
	free(T);
}

/*
DM2 using binomial distribution to calculate the P-value
the percentile is chose to be upper 50%
*/

void
DeviationMeasure2(int n)
{
    int       i, ii, j, d, N, L, m, N_;
    double    p_value, T_, nu;
    double fenweidian = 533.5;
    double    pi = 0.4978132181;

    BitSequence  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
    const int M = 500;
    N = (int)floor(n / M);
    if (((B_ = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
        ((C = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
        ((P = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL) ||
        ((T = (BitSequence *)calloc(M + 1, sizeof(BitSequence))) == NULL)) {
        printf("Insufficient Memory for Work Space:: Deviation Measure II Test\n");
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

    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tDEVIATION MEASURE II\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tNote: %d bits were discarded!\n", n%M);
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "  Frequency     CHI2      P-value\n");

    nu = 0.00;
    for (ii = 0; ii<N; ii++) {
        T_ = 0.0;
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
                    L = N_ + 1 - L;
                    m = N_;
                    for (i = 0; i<M + 1; i++)
                        B_[i] = T[i];
                }
            }
            T_ += (L - (N_ + 1) / 2.0)*(L - (N_ + 1) / 2.0);
            N_++;
        }



        if (T_ > fenweidian)
            nu += 1.0;
    }
    double p = nu / N;
    double z = (p - pi) / sqrt(pi * (1 - pi) / N);

    fprintf(stats[TEST_DEVIATIONMEARSURE2], "    %4d    ", (int)nu);
    p_value = cephes_erfc(fabs(z) / sqrt(2));

    fprintf(stats[TEST_DEVIATIONMEARSURE2], "%9.6f    %9.6f\n", z, p_value); fflush(stats[TEST_DEVIATIONMEARSURE2]);
    fprintf(results[TEST_DEVIATIONMEARSURE2], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE2]);

    free(B_);
    free(P);
    free(C);
    free(T);
}

void
DeviationMeasure1CS1(int n)
{
	int       i, ii, N, K = 6;
	double    p_value, nu[7], chi2;
    double fenweidian[6] = { 335.0,
                             351.0,
                             364.0,
                             377.0,
                             393.0,
                             415.0 };
    double    pi[7]={ 0.1400018365,
                      0.1422095933,
                      0.1394350127,
                      0.1403124093,
                      0.1503319371,
                      0.1444503240,
                      0.1432588871 };
	
	double 	  *counter;
	unsigned int    power2[32], lenS, bitLenS;

	const int M = 500;
	N = (int)floor(n/M); //分块块数
	lenS = M; 
	bitLenS = (lenS + 31) / 32;
	power2[0] = 1;
	for (int i = 1; i < 32; i++)
		power2[i] = 2 * power2[i - 1];
	if( (counter = (double*)malloc(sizeof(double) * N)) == NULL ){
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		return;		
	}

	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tDEVIATION MEASURE I CS\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tNote: %d bits were discarded!\n", n%M);

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
		double T_ ;
		// int word, bitPos;
		n = lenC = lenD = 0;
		m = -1;
		T_ = 0.0;
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
			T_ += fabs(lenC - (n + 1)/2.0);
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
		double T_=counter[ii];
		i = 0;
        while (i<6 && T_ >= fenweidian[i])
            i++;
        nu[i] += 1.0;
	}
	
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ )
		fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "%4d ", (int)nu[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_DEVIATIONMEARSURE1CS]);
	fprintf(results[TEST_DEVIATIONMEARSURE1CS], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE1CS]);
}

void
DeviationMeasure1_1(int n)
{
	int       ii, N;
	double    p_value, nu;
    double    fenweidian = 371.0;
    double    pi = 0.4909290857;
	
	double 	  *counter;
	unsigned int    power2[32], lenS, bitLenS;

	const int M = 500;
	N = (int)floor(n/M); //分块块数
	lenS = M; 
	bitLenS = (lenS + 31) / 32;
	power2[0] = 1;
	for (int i = 1; i < 32; i++)
		power2[i] = 2 * power2[i - 1];
	if( (counter = (double*)malloc(sizeof(double) * N)) == NULL ){
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		return;		
	}

    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tDEVIATION MEASURE I\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tNote: %d bits were discarded!\n", n%M);
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "  Frequency     CHI2      P-value\n");

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
		double T_ ;
		// int word, bitPos;
		n = lenC = lenD = 0;
		m = -1;
		T_ = 0.0;
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
			T_ += fabs(lenC - (n + 1) / 2.0);
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

	nu=0.00;
	for( ii=0; ii<N; ii++){
		double T_=counter[ii];
        if(T_ > fenweidian)
            nu += 1.0;
	}

    double p = nu / N;
    double z = (p - pi) / sqrt(pi * (1 - pi) / N);
    
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "    %4d    ", (int)nu);
    p_value = cephes_erfc(fabs(z) / sqrt(2));

    fprintf(stats[TEST_DEVIATIONMEARSURE1], "%9.6f    %9.6f\n", z, p_value); fflush(stats[TEST_DEVIATIONMEARSURE1]);
    fprintf(results[TEST_DEVIATIONMEARSURE1], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE1]);
}

void
DeviationMeasure2CS1(int n)
{
	int       i, ii, N, K = 6;
	double    p_value, nu[7], chi2;
    double fenweidian[6] = {
        422.5,
        469.5,
        511.5,
        556.5,
        612.5,
        699.5
    };
	double    pi[7]={  
		0.1404715615,
        0.1427569687,
        0.1430949506,
        0.1436457382,
        0.1440210422,
        0.1423812606,
        0.1436284782 };
	
	double 	  *counter;
	unsigned int    power2[32], lenS, bitLenS;

	const int M = 500;
	N = (int)floor(n/M); //分块块数
	lenS = M; 
	bitLenS = (lenS + 31) / 32;
	power2[0] = 1;
	for (int i = 1; i < 32; i++)
		power2[i] = 2 * power2[i - 1];
	if( (counter = (double*)malloc(sizeof(double) * N)) == NULL ){
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		return;		
	}

	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tDEVIATION MEASURE II CS\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tNote: %d bits were discarded!\n", n%M);

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
		double T_ ;
		// int word, bitPos;
		n = lenC = lenD = 0;
		m = -1;
		T_ = 0.0;
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
			T_ += (lenC - (n+1)/2.0)*(lenC - (n+1)/2.0);
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
		double T_=counter[ii];
		i = 0;
        while (i<6 && T_ >= fenweidian[i])
            i++;
        nu[i] += 1.0;
	}
	
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ )
		fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "%4d ", (int)nu[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_DEVIATIONMEARSURE2CS]);
	fprintf(results[TEST_DEVIATIONMEARSURE2CS], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE2CS]);
}

void
DeviationMeasure2_1(int n)
{
	int       ii, N;
	double    p_value, nu;
    double    fenweidian = 533.5;
    double    pi = 0.4978132181;
	
	double 	  *counter;
	unsigned int    power2[32], lenS, bitLenS;

	const int M = 500;
	N = (int)floor(n/M); //分块块数
	lenS = M; 
	bitLenS = (lenS + 31) / 32;
	power2[0] = 1;
	for (int i = 1; i < 32; i++)
		power2[i] = 2 * power2[i - 1];
	if( (counter = (double*)malloc(sizeof(double) * N)) == NULL ){
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		return;		
	}

    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tDEVIATION MEASURE II\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tNote: %d bits were discarded!\n", n%M);
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "  Frequency     CHI2      P-value\n");

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
		double T_ ;
		// int word, bitPos;
		n = lenC = lenD = 0;
		m = -1;
		T_ = 0.0;
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
			T_ += (lenC - (n + 1) / 2.0)*(lenC - (n + 1) / 2.0);
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

	nu=0.00;
	for( ii=0; ii<N; ii++){
		double T_=counter[ii];
        if(T_ > fenweidian)
            nu += 1.0;
	}
	
    double p = nu / N;
    double z = (p - pi) / sqrt(pi * (1 - pi) / N);

    fprintf(stats[TEST_DEVIATIONMEARSURE2], "    %4d    ", (int)nu);
    p_value = cephes_erfc(fabs(z) / sqrt(2));

    fprintf(stats[TEST_DEVIATIONMEARSURE2], "%9.6f    %9.6f\n", z, p_value); fflush(stats[TEST_DEVIATIONMEARSURE2]);
    fprintf(results[TEST_DEVIATIONMEARSURE2], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE2]);
}

void
DeviationMeasureTotal1(int n)
{
	int       i, ii, N, K = 6;
	double    p_value, chi2, p, z;
	double    nucs[7];
	double    nu;
    
	double    fenweidian_cs1[6] = { 335.0, 351.0, 364.0, 377.0, 393.0, 415.0 };
    double    pi_cs1[7]={ 0.1400018365, 0.1422095933, 0.1394350127, 0.1403124093,
							 0.1503319371, 0.1444503240, 0.1432588871 };
    double    fenweidian_cs2[6] = { 422.5, 469.5, 511.5, 556.5, 612.5, 699.5};
	double    pi_cs2[7]={ 0.1404715615, 0.1427569687, 0.1430949506, 0.1436457382,
        					0.1440210422, 0.1423812606, 0.1436284782 };
    double    fenweidian_1 = 371.0;
    double    pi_1 = 0.4909290857;
    double    fenweidian_2 = 533.5;
    double    pi_2 = 0.4978132181;

	double 	  *counter[2];
	unsigned int    power2[32], lenS, bitLenS;
	const int M = 500;
	N = (int)floor(n/M); //分块块数
	lenS = M; 
	bitLenS = (lenS + 31) / 32;
	power2[0] = 1;
	for (int i = 1; i < 32; i++)
		power2[i] = 2 * power2[i - 1];
    for( i = 0; i < 2; i++){
        if( (counter[i] = (double*)malloc(sizeof(double) * N)) == NULL ){
            printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
            return;		
        }        
    }

	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tDEVIATION MEASURE I CS\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "\tNote: %d bits were discarded!\n", n%M);

    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tDEVIATION MEASURE I\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "\tNote: %d bits were discarded!\n", n%M);
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "  Frequency     CHI2      P-value\n");

	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tDEVIATION MEASURE II CS\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "  C0   C1   C2   C3   C4   C5   C6     CHI2    P-value\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "-----------------------------------------------------\n");
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "\tNote: %d bits were discarded!\n", n%M);

    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tDEVIATION MEASURE II\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tM (substring length)     = %d\n", M);
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tN (number of substrings) = %d\n", N);
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "\tNote: %d bits were discarded!\n", n%M);
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "-----------------------------------------------------\n");
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "  Frequency     CHI2      P-value\n");

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
		double T_1, T_2 ;
		// int word, bitPos;
		n = lenC = lenD = 0;
		m = -1;
		T_1 = 0.0;
		T_2 = 0.0;
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
			T_1 += fabs(lenC - (n + 1)/2.0);
			T_2 += (lenC - (n+1)/2.0)*(lenC - (n+1)/2.0);
			n++;
		}

		//存入
		counter[0][ii] = T_1;
		counter[1][ii] = T_2;

		free(bitC);
		free(--bitD);
		free(bitTmp);

		for (int i = 0; i < 32; i++)
			free(bitRS[i]);
	}


	for ( i=0; i<K+1; i++ )
		nucs[i] = 0.00;
	for( ii=0; ii<N; ii++){
		double T_=counter[0][ii];
		i = 0;
        while (i<6 && T_ >= fenweidian_cs1[i])
            i++;
        nucs[i] += 1.0;
	}
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ )
		fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "%4d ", (int)nucs[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nucs[i]-N*pi_cs1[i], 2) / (N*pi_cs1[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);
	fprintf(stats[TEST_DEVIATIONMEARSURE1CS], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_DEVIATIONMEARSURE1CS]);
	fprintf(results[TEST_DEVIATIONMEARSURE1CS], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE1CS]);


	nu=0.00;
	for( ii=0; ii<N; ii++){
		double T_=counter[0][ii];
        if(T_ > fenweidian_1)
            nu += 1.0;
	}
    p = nu / N;
    z = (p - pi_1) / sqrt(pi_1 * (1 - pi_1) / N);
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "    %4d    ", (int)nu);
    p_value = cephes_erfc(fabs(z) / sqrt(2));
    fprintf(stats[TEST_DEVIATIONMEARSURE1], "%9.6f    %9.6f\n", z, p_value); fflush(stats[TEST_DEVIATIONMEARSURE1]);
    fprintf(results[TEST_DEVIATIONMEARSURE1], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE1]);


	for ( i=0; i<K+1; i++ )
		nucs[i] = 0.00;
	for( ii=0; ii<N; ii++){
		double T_=counter[1][ii];
		i = 0;
        while (i<6 && T_ >= fenweidian_cs2[i])
            i++;
        nucs[i] += 1.0;
	}
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ )
		fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "%4d ", (int)nucs[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nucs[i]-N*pi_cs2[i], 2) / (N*pi_cs2[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);
	fprintf(stats[TEST_DEVIATIONMEARSURE2CS], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_DEVIATIONMEARSURE2CS]);
	fprintf(results[TEST_DEVIATIONMEARSURE2CS], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE2CS]);


	nu=0.00;
	for( ii=0; ii<N; ii++){
		double T_=counter[1][ii];
        if(T_ > fenweidian_2)
            nu += 1.0;
	}
    p = nu / N;
    z = (p - pi_2) / sqrt(pi_2 * (1 - pi_2) / N);
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "    %4d    ", (int)nu);
    p_value = cephes_erfc(fabs(z) / sqrt(2));
    fprintf(stats[TEST_DEVIATIONMEARSURE2], "%9.6f    %9.6f\n", z, p_value); fflush(stats[TEST_DEVIATIONMEARSURE2]);
    fprintf(results[TEST_DEVIATIONMEARSURE2], "%f\n", p_value); fflush(results[TEST_DEVIATIONMEARSURE2]);
}
