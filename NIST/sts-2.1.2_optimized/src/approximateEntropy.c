#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <arpa/inet.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                A P P R O X I M A T E  E N T R O P Y   T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
ApproximateEntropy(int m, int n)
{
	int				i, j, k, r, blockSize, seqLength, powLen, index;
	double			sum, numOfBlocks, ApEn[2], apen, chi_squared, p_value;
	unsigned int	*P;
	
	fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);

	seqLength = n;
	r = 0;
	
	for ( blockSize=m; blockSize<=m+1; blockSize++ ) {
		if ( blockSize == 0 ) {
			ApEn[0] = 0.00;
			r++;
		}
		else {
			numOfBlocks = (double)seqLength;
			powLen = (int)pow(2, blockSize+1)-1;
			if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== NULL ) {
				fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
				return;
			}
			for ( i=1; i<powLen-1; i++ )
				P[i] = 0;
			for ( i=0; i<numOfBlocks; i++ ) { /* COMPUTE FREQUENCY */
				k = 1;
				for ( j=0; j<blockSize; j++ ) {
					k <<= 1;
					if ( (int)epsilon[(i+j) % seqLength] == 1 )
						k++;
				}
				P[k-1]++;
			}
			/* DISPLAY FREQUENCY */
			sum = 0.0;
			index = (int)pow(2, blockSize)-1;
			for ( i=0; i<(int)pow(2, blockSize); i++ ) {
				if ( P[index] > 0 )
					sum += P[index]*log(P[index]/numOfBlocks);
				index++;
			}
			sum /= numOfBlocks;
			ApEn[r] = sum;
			r++;
			free(P);
		}
	}
	apen = ApEn[0] - ApEn[1];
	
	chi_squared = 2.0*seqLength*(log(2) - apen);
	p_value = cephes_igamc(pow(2, m-1), chi_squared/2.0);
	
	fprintf(stats[TEST_APEN], "\t\t(b) n (sequence length) = %d\n", seqLength);
	fprintf(stats[TEST_APEN], "\t\t(c) Chi^2               = %f\n", chi_squared);
	fprintf(stats[TEST_APEN], "\t\t(d) Phi(m)	       = %f\n", ApEn[0]);
	fprintf(stats[TEST_APEN], "\t\t(e) Phi(m+1)	       = %f\n", ApEn[1]);
	fprintf(stats[TEST_APEN], "\t\t(f) ApEn                = %f\n", apen);
	fprintf(stats[TEST_APEN], "\t\t(g) Log(2)              = %f\n", log(2.0));
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");

	if ( m > (int)(log(seqLength)/log(2)-5) ) {
		fprintf(stats[TEST_APEN], "\t\tNote: The blockSize = %d exceeds recommended value of %d\n", m,
			MAX(1, (int)(log(seqLength)/log(2)-5)));
		fprintf(stats[TEST_APEN], "\t\tResults are inaccurate!\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	}
	
	fprintf(stats[TEST_APEN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_APEN]);
	fprintf(results[TEST_APEN], "%f\n", p_value); fflush(results[TEST_APEN]);
}

// 近似熵检测
void
ApproximateEntropy1(int m, int n)
{
	int				i, j, k, r, blockSize, seqLength, powLen, index, byteNum, M;
	double			sum, numOfBlocks, ApEn[2], apen, chi_squared, p_value;
	unsigned int	*P;
	
	fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);

	seqLength = n;
	r = 0;
	
	for ( blockSize=m; blockSize<=m+1; blockSize++ ) {
		if ( blockSize == 0 ) {
			ApEn[0] = 0.00;
			r++;
		}
		else {
			numOfBlocks = (double)seqLength;
			
			if(blockSize <= 32){
				byteNum = numOfBlocks / 8; 
				powLen = (int)pow(2, blockSize);
				M = powLen - 1;
				// 数组P存储所有可能重叠的2^m的模式的频数数组
				if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== NULL ) {
					fprintf(stats[TEST_SERIAL], "Serial Test:  Insufficient memory available.\n");
					fflush(stats[TEST_SERIAL]);
					return ;
				}
				for ( i=0; i<powLen; i++ )
					P[i] = 0;	  /* INITIALIZE NODES */

				// int *container = (int*)calloc(n,sizeof(int));
				// #pragma omp parallel for 

				for ( int i=0; i<byteNum; ++i ) {		 /* COMPUTE FREQUENCY */
					for (int j = 0; j < 8; ++j) {
						int k = 8 - j - blockSize;
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
							temp = ((htonl(((int*)(byteEpsilon+i))[0]) << (-24 - k))|(byteEpsilon[i+4]>>(32+k))) & M;
						}
						P[temp] += 1;
					}
				}

				//处理尾巴
				if (tp.n % 8 != 0) {
					int i = byteNum;
					for (int j = 0; j < tp.n - 8 * i; j++) {
						int k = 8 - j - blockSize;
						int temp;
						if (k >= 0) {
							// 涉及1byte
							temp = (byteEpsilon[i] >> k) & M;
						}
						else if (k >= -8) {
							//涉及2bytes
							temp = ((byteEpsilon[i] << (-k)) | (byteEpsilon[i + 1] >> (8 + k))) & M;
						}
						else if (k >= -24) {
							//涉及3bytes
							temp = (htonl(((int*)(byteEpsilon+i))[0]) >> (24 + k)) & M;
						}
						else if (k >= -31) {
							//涉及4bytes
							temp = ((htonl(((int*)(byteEpsilon+i))[0]) << (-24 - k)) | (byteEpsilon[i + 4] >> (32 + k))) & M;
						}
						P[temp] += 1;
					}
				}
				
				/* DISPLAY FREQUENCY */
				sum = 0.0;
				for ( i=0; i<powLen; i++ ) {
					if ( P[i] > 0 )
						sum += P[i]*log(P[i]/numOfBlocks);
				}	
			}else{
				powLen = (int)pow(2, blockSize+1)-1;
				if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== NULL ) {
					fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
					return;
				}
				for ( i=1; i<powLen-1; i++ )
					P[i] = 0;
				for ( i=0; i<numOfBlocks; i++ ) { /* COMPUTE FREQUENCY */
					k = 1;
					for ( j=0; j<blockSize; j++ ) {
						k <<= 1;
						if ( (int)epsilon[(i+j) % seqLength] == 1 )
							k++;
					}
					P[k-1]++;
				}
				/* DISPLAY FREQUENCY */
				sum = 0.0;
				index = (int)pow(2, blockSize)-1;
				for ( i=0; i<(int)pow(2, blockSize); i++ ) {
					if ( P[index] > 0 )
						sum += P[index]*log(P[index]/numOfBlocks);
					index++;
				}
			}
			
			sum /= numOfBlocks;
			ApEn[r] = sum;
			r++;
			free(P);
		}
	}
	apen = ApEn[0] - ApEn[1];
	
	chi_squared = 2.0*seqLength*(log(2) - apen);
	p_value = cephes_igamc(pow(2, m-1), chi_squared/2.0);
	
	fprintf(stats[TEST_APEN], "\t\t(b) n (sequence length) = %d\n", seqLength);
	fprintf(stats[TEST_APEN], "\t\t(c) Chi^2               = %f\n", chi_squared);
	fprintf(stats[TEST_APEN], "\t\t(d) Phi(m)	       = %f\n", ApEn[0]);
	fprintf(stats[TEST_APEN], "\t\t(e) Phi(m+1)	       = %f\n", ApEn[1]);
	fprintf(stats[TEST_APEN], "\t\t(f) ApEn                = %f\n", apen);
	fprintf(stats[TEST_APEN], "\t\t(g) Log(2)              = %f\n", log(2.0));
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");

	if ( m > (int)(log(seqLength)/log(2)-5) ) {
		fprintf(stats[TEST_APEN], "\t\tNote: The blockSize = %d exceeds recommended value of %d\n", m,
			MAX(1, (int)(log(seqLength)/log(2)-5)));
		fprintf(stats[TEST_APEN], "\t\tResults are inaccurate!\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	}
	
	fprintf(stats[TEST_APEN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_APEN]);
	fprintf(results[TEST_APEN], "%f\n", p_value); fflush(results[TEST_APEN]);
}