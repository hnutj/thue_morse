#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <arpa/inet.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
               O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double	Pr(int u, double eta);

void
OverlappingTemplateMatchings(int m, int n)
{
	int				i, k, match;
	double			W_obs, eta, sum, chi2, p_value, lambda;
	int				M, N, j, K = 5;
	unsigned int	nu[6] = { 0, 0, 0, 0, 0, 0 };
	//double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 };
	double			pi[6] = { 0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865 };
	BitSequence		*sequence;

	M = 1032;
	N = n/M;
	
	if ( (sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == NULL ) {
		fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\tTEMPLATE DEFINITION:  Insufficient memory, Overlapping Template Matchings test aborted!\n");
	}
	else
		for ( i=0; i<m; i++ )
			sequence[i] = 1;
	
	lambda = (double)(M-m+1)/pow(2,m);
	eta = lambda/2.0;
	sum = 0.0;
	for ( i=0; i<K; i++ ) {			/* Compute Probabilities */
		pi[i] = Pr(i, eta);
		sum += pi[i];
	}
	pi[K] = 1 - sum;

	for ( i=0; i<N; i++ ) {
		W_obs = 0;
		for ( j=0; j<M-m+1; j++ ) {
			match = 1;
			for ( k=0; k<m; k++ ) {
				if ( sequence[k] != epsilon[i*M+j+k] )
					match = 0;
			}
			if ( match == 1 )
				W_obs++;
		}
		if ( W_obs <= 4 )
			nu[(int)W_obs]++;
		else
			nu[K]++;
	}
	sum = 0;
	chi2 = 0.0;                                   /* Compute Chi Square */
	for ( i=0; i<K+1; i++ ) {
		chi2 += pow((double)nu[i] - (double)N*pi[i], 2)/((double)N*pi[i]);
		sum += nu[i];
	}
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t(a) n (sequence_length)      = %d\n", n);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(b) m (block length of 1s)   = %d\n", m);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(c) M (length of substring)  = %d\n", M);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(d) N (number of substrings) = %d\n", N);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(f) eta                      = %f\n", eta);
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t   F R E Q U E N C Y\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t  0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t%3d %3d %3d %3d %3d %3d  %f ",
		nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_OVERLAPPING], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

	free(sequence);
	fprintf(stats[TEST_OVERLAPPING], "%f %s\n\n", p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS"); fflush(stats[TEST_OVERLAPPING]);
	fprintf(results[TEST_OVERLAPPING], "%f\n", p_value); fflush(results[TEST_OVERLAPPING]);
}

void
OverlappingTemplateMatchings1(int m, int n)
{
	// double          W_obs;
	// int             i, j, k, match;
	int				i;
	double			eta, sum, chi2, p_value, lambda;
	int				M, N, K = 5;
	unsigned int	nu[6] = { 0, 0, 0, 0, 0, 0 };
	//double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 };
	double			pi[6] = { 0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865 };
	BitSequence		*sequence;

	M = 1032;
	N = n/M;
	
	if ( (sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == NULL ) {
		fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\tTEMPLATE DEFINITION:  Insufficient memory, Overlapping Template Matchings test aborted!\n");
	}
	else
		for ( i=0; i<m; i++ )
			sequence[i] = 1;
	
	lambda = (double)(M-m+1)/pow(2,m);
	eta = lambda/2.0;
	sum = 0.0;
	for ( i=0; i<K; i++ ) {			/* Compute Probabilities */
		pi[i] = Pr(i, eta);
		sum += pi[i];
	}
	pi[K] = 1 - sum;

	unsigned short mask = pow(2, m) - 1;
	unsigned short templates = mask;
	int *Wj = (int*)calloc(N,sizeof(int));
	
	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		//每个子块的索引区域：j*M~(j+1)*M-m
		int beginPos = i * M, endPos = beginPos + M - m;
		int beginByte = beginPos / 8, beginOff = beginPos % 8;
		int endByte = endPos / 8, endOff = endPos % 8;

		int W_obs = 0;
		for (int p = beginByte; p <= endByte; p++) {
			int pos1 = (p == beginByte ? beginOff : 0);
			int pos2 = (p == endByte ? endOff : 7);
			for (int t = pos1; t <= pos2; t++) {
				unsigned short test;
				int k = 8 - t - m;
				if (k >= 0) {
					// 涉及1byte
					test = (byteEpsilon[p] >> k) & mask;
				}
				else if (k >= -8) {
					//涉及2bytes
					test = ((byteEpsilon[p] << (-k)) | (byteEpsilon[p + 1] >> (8 + k))) & mask;
				}
				else if (k >= -15) {
					//涉及3bytes
					test = (htonl(((int*)(byteEpsilon + p))[0]) >> (24 + k)) & mask;
				}
				if (test == templates) {
					W_obs++;
				}
			}	
		}
		Wj[i]=W_obs;
	}
	for(i=0;i<N;i++){
		int W_obs = Wj[i];
		if ( W_obs <= 4 )
			nu[W_obs]++;
		else
			nu[K]++;		
	}

	sum = 0;
	chi2 = 0.0;                                   /* Compute Chi Square */
	for ( i=0; i<K+1; i++ ) {
		chi2 += pow((double)nu[i] - (double)N*pi[i], 2)/((double)N*pi[i]);
		sum += nu[i];
	}
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t(a) n (sequence_length)      = %d\n", n);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(b) m (block length of 1s)   = %d\n", m);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(c) M (length of substring)  = %d\n", M);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(d) N (number of substrings) = %d\n", N);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(f) eta                      = %f\n", eta);
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t   F R E Q U E N C Y\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t  0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t%3d %3d %3d %3d %3d %3d  %f ",
		nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_OVERLAPPING], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

	free(sequence);
	fprintf(stats[TEST_OVERLAPPING], "%f %s\n\n", p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS"); fflush(stats[TEST_OVERLAPPING]);
	fprintf(results[TEST_OVERLAPPING], "%f\n", p_value); fflush(results[TEST_OVERLAPPING]);
}

double
Pr(int u, double eta)
{
	int		l;
	double	sum, p;
	
	if ( u == 0 )
		p = exp(-eta);
	else {
		sum = 0.0;
		for ( l=1; l<=u; l++ )
			sum += exp(-eta-u*log(2)+l*log(eta)-cephes_lgam(l+1)+cephes_lgam(u)-cephes_lgam(l)-cephes_lgam(u-l+1));
		p = sum;
	}
	return p;
}
