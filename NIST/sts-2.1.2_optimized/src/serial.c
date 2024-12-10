#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <arpa/inet.h>
#include "../include/externs.h"
#include "../include/cephes.h"  

double psi2(int m, int n);
double myPsi2(int m, int n);


void
Serial(int m, int n)
{
	double	p_value1, p_value2, psim0, psim1, psim2, del1, del2;
	
	psim0 = psi2(m, n);
	psim1 = psi2(m-1, n);
	psim2 = psi2(m-2, n);
	
	del1 = psim0 - psim1;
	del2 = psim0 - 2.0*psim1 + psim2;
	p_value1 = cephes_igamc(pow(2, m-1)/2, del1/2.0);
	p_value2 = cephes_igamc(pow(2, m-2)/2, del2/2.0);
	
	fprintf(stats[TEST_SERIAL], "\t\t\t       SERIAL TEST\n");
	fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_SERIAL], "\t\t COMPUTATIONAL INFORMATION:		  \n");
	fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_SERIAL], "\t\t(a) Block length    (m) = %d\n", m);
	fprintf(stats[TEST_SERIAL], "\t\t(b) Sequence length (n) = %d\n", n);
	fprintf(stats[TEST_SERIAL], "\t\t(c) Psi_m               = %f\n", psim0);
	fprintf(stats[TEST_SERIAL], "\t\t(d) Psi_m-1             = %f\n", psim1);
	fprintf(stats[TEST_SERIAL], "\t\t(e) Psi_m-2             = %f\n", psim2);
	fprintf(stats[TEST_SERIAL], "\t\t(f) Del_1               = %f\n", del1);
	fprintf(stats[TEST_SERIAL], "\t\t(g) Del_2               = %f\n", del2);
	fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_SERIAL], "%s\t\tp_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
	fprintf(results[TEST_SERIAL], "%f\n", p_value1);

	fprintf(stats[TEST_SERIAL], "%s\t\tp_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2); fflush(stats[TEST_SERIAL]);
	fprintf(results[TEST_SERIAL], "%f\n", p_value2); fflush(results[TEST_SERIAL]);
}

// 序列检验
void
Serial1(int m, int n)
{
	double	p_value1, p_value2, psim0, psim1, psim2, del1, del2;
	
	// psim0 = psi2(m, n);
	// psim1 = psi2(m-1, n);
	// psim2 = psi2(m-2, n);

	psim0 = myPsi2(m, n);
	psim1 = myPsi2(m-1, n);
	psim2 = myPsi2(m-2, n);

	del1 = psim0 - psim1;
	del2 = psim0 - 2.0*psim1 + psim2;
	p_value1 = cephes_igamc(pow(2, m-1)/2, del1/2.0);
	p_value2 = cephes_igamc(pow(2, m-2)/2, del2/2.0);
	
	fprintf(stats[TEST_SERIAL], "\t\t\t       SERIAL TEST\n");
	fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_SERIAL], "\t\t COMPUTATIONAL INFORMATION:		  \n");
	fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_SERIAL], "\t\t(a) Block length    (m) = %d\n", m);
	fprintf(stats[TEST_SERIAL], "\t\t(b) Sequence length (n) = %d\n", n);
	fprintf(stats[TEST_SERIAL], "\t\t(c) Psi_m               = %f\n", psim0);
	fprintf(stats[TEST_SERIAL], "\t\t(d) Psi_m-1             = %f\n", psim1);
	fprintf(stats[TEST_SERIAL], "\t\t(e) Psi_m-2             = %f\n", psim2);
	fprintf(stats[TEST_SERIAL], "\t\t(f) Del_1               = %f\n", del1);
	fprintf(stats[TEST_SERIAL], "\t\t(g) Del_2               = %f\n", del2);
	fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_SERIAL], "%s\t\tp_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
	fprintf(results[TEST_SERIAL], "%f\n", p_value1);

	fprintf(stats[TEST_SERIAL], "%s\t\tp_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2); fflush(stats[TEST_SERIAL]);
	fprintf(results[TEST_SERIAL], "%f\n", p_value2); fflush(results[TEST_SERIAL]);
}

// n是序列总长度，m是子序列的长度（块长）
double
psi2(int m, int n)
{
	int				i, j, k, powLen;
	double			sum, numOfBlocks;
	unsigned int	*P;
	
	if ( (m == 0) || (m == -1) )
		return 0.0;
	numOfBlocks = n;
	powLen = (int)pow(2, m+1)-1;
	// 数组P存储所有可能重叠的2^m的模式的频数数组
	if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== NULL ) {
		fprintf(stats[TEST_SERIAL], "Serial Test:  Insufficient memory available.\n");
		fflush(stats[TEST_SERIAL]);
		return 0.0;
	}
	for ( i=1; i<powLen-1; i++ )
		P[i] = 0;	  /* INITIALIZE NODES */
	for ( i=0; i<numOfBlocks; i++ ) {		 /* COMPUTE FREQUENCY */
		k = 1;
		for ( j=0; j<m; j++ ) {
			if ( epsilon[(i+j)%n] == 0 )
				k *= 2;
			else if ( epsilon[(i+j)%n] == 1 )
				k = 2*k+1;
		}
		P[k-1]++;
	}
	sum = 0.0;
	for ( i=(int)pow(2, m)-1; i<(int)pow(2, m+1)-1; i++ )
		sum += pow(P[i], 2);
	sum = (sum * pow(2, m)/(double)n) - (double)n;
	free(P);
	
	return sum;
}

// n是序列总长度，m是子序列的长度,要求m<=8
double myPsi2(int m, int n){
	int				i, M, powLen, byteNum;
	double			sum;
	unsigned int	*P;
	
	if ( (m == 0) || (m == -1) )
		return 0.0;
	byteNum = n / 8; //这边以字节为单位
	powLen = (int)pow(2, m);
	M = powLen - 1;
	// 数组P存储所有可能重叠的2^m的模式的频数数组
	if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== NULL ) {
		fprintf(stats[TEST_SERIAL], "Serial Test:  Insufficient memory available.\n");
		fflush(stats[TEST_SERIAL]);
		return 0.0;
	}
	for ( i=0; i<powLen; i++ )
		P[i] = 0;	  /* INITIALIZE NODES */

	// int *container = (int*)calloc(n,sizeof(int));
	// #pragma omp parallel for 
	
	for ( int i=0; i<byteNum; ++i ) {		 /* COMPUTE FREQUENCY */
		for ( int j = 0; j < 8; ++j) {
			int k = 8 - j - m;
			int temp;
			if (k >= 0) {
				// 涉及1byte
				temp = (byteEpsilon[i] >> k) & M;				
			}
			else if(k >= -8){
				//涉及2bytes
				temp = ((byteEpsilon[i] << (-k)) | (byteEpsilon[i + 1] >> (8 + k))) & M;
			}
			else if (k >= -24) {
				//涉及3bytes
				temp = (htonl(((int*)(byteEpsilon+i))[0]) >> (24 + k)) & M;
			}
			else if (k >= -31) {
				//涉及4bytes
				temp = ((htonl(((int*)(byteEpsilon+i))[0]) << (-24 - k))|(byteEpsilon[i + 4] >> (32+k))) & M;
			}
			P[temp] += 1;
		}
	}

	//处理尾巴
	if (tp.n % 8 != 0) {
		int i = byteNum;
		for (int j = 0; j < tp.n - 8 * i; j++) {
			int k = 8 - j - m;
			int temp;
			if (k >= 0) {
				// 涉及1byte
				temp = (byteEpsilon[i] >> k) & M;				
			}
			else if(k >= -8){
				//涉及2bytes
				temp = ((byteEpsilon[i] << (-k)) | (byteEpsilon[i + 1] >> (8 + k))) & M;
			}
			else if (k >= -24) {
				//涉及3bytes
				temp = (htonl(((int*)(byteEpsilon+i))[0]) >> (24 + k)) & M;
			}
			else if (k >= -31) {
				//涉及4bytes
				temp = ((htonl(((int*)(byteEpsilon+i))[0]) << (-24 - k))|(byteEpsilon[i + 4] >> (32+k))) & M;
			}
			P[temp] += 1;
		}
	}

	sum = 0.0;
	for ( i=0; i<powLen; i++ )
		sum += pow(P[i], 2);
	sum = (sum * pow(2, m)/(double)n) - (double)n;
	free(P);
	
	return sum;
}