#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <arpa/inet.h>
#include "../include/cephes.h"
#include "../include/externs.h"

void
PokerDetect(int m, int n)
{
	// for(int i=0;i<n;i++){
    //     printf("%1d ",epsilon[i]);
    // }
    // printf("\n");
 	int				i, j, k, N, powLen;
	double			sum, V, p_value;
	unsigned int	*P;
    if ( (m == 0) || (m == -1) )
		return;
    N = n / m;
	powLen = (int)pow(2, m+1)-1;
	// 数组P存储所有可能重叠的2^m的模式的频数数组
	if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== NULL ) {
		fprintf(stats[TEST_POKERDETECT], "PokerDetect Test:  Insufficient memory available.\n");
		fflush(stats[TEST_POKERDETECT]);
		return;
	}

	for ( i=1; i<powLen-1; i++ )
		P[i] = 0;	  /* INITIALIZE NODES */
    for ( i=0; i<N*m; i+=m ){
        k = 1;
		for ( j=0; j<m; j++ ) {
			if ( epsilon[(i+j)%n] == 0 )
				k *= 2;
			else if ( epsilon[(i+j)%n] == 1 )
				k = 2*k+1;
		}
		P[k-1]++;
		// printf("i=%d, k=%d\n",i,k);
    }

    sum = 0.0;
	for ( i=(int)pow(2, m)-1; i<powLen; i++ )
		sum += pow(P[i], 2);
	// printf("sum=%f\n",sum);
    V = (pow(2, m) / N) * sum - N;
    p_value = cephes_igamc((pow(2, m)-1) / 2, V / 2);

    fprintf(stats[TEST_POKERDETECT], "\t\t\t       POKERDETECT TEST\n");
	fprintf(stats[TEST_POKERDETECT], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_POKERDETECT], "\t\t COMPUTATIONAL INFORMATION:		  \n");
	fprintf(stats[TEST_POKERDETECT], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_POKERDETECT], "\t\t(a) Block length    (m) = %d\n", m);
	fprintf(stats[TEST_POKERDETECT], "\t\t(b) Sequence length (n) = %d\n", n);
	fprintf(stats[TEST_POKERDETECT], "\t\t(c) V                   = %f\n", V);
	fprintf(stats[TEST_POKERDETECT], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_POKERDETECT], "%s\t\tp_value = %f\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_POKERDETECT], "%f\n", p_value);
}

void
PokerDetect1(int m, int n)
{
	// for(int i=0;i<n;i++){
    //     printf("%1d ",epsilon[i]);
    // }
    // printf("\n");
 	int				i, j, k, N, powLen, M, temp;
	double			sum, V, p_value;
	unsigned int	*P;
    if ( (m == 0) || (m == -1) )
		return;
    N = n / m;
	powLen = (int)pow(2, m);
	M = powLen - 1;
	// 数组P存储所有可能重叠的2^m的模式的频数数组
	if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== NULL ) {
		fprintf(stats[TEST_POKERDETECT], "PokerDetect Test:  Insufficient memory available.\n");
		fflush(stats[TEST_POKERDETECT]);
		return;
	}

	for ( i = 0; i < powLen; i++ ) {
		P[i] = 0;
	}

	if(m == 8){
		int byteNum = n / 8;
		for( i = 0; i < byteNum; i++){
			P[byteEpsilon[i]]++;
		}
	}
	else{
		for (int pos = 0; pos < N * m; pos += m) {
			i = pos / 8;
			j = pos % 8;
			k = 8 - j - m;
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
				temp = (htonl(((int*)(byteEpsilon + i))[0]) >> (24 + k)) & M;
			}
			else if (k >= -31) {
				//涉及4bytes
				temp = ((htonl(((int*)(byteEpsilon + i))[0]) << (-24 - k)) | (byteEpsilon[i + 4] >> (32 + k))) & M;
			}
			P[temp] += 1;
		}		
	}


    sum = 0.0;
	for ( i = 0; i < powLen; i++ )
		sum += pow(P[i], 2);
	// printf("sum=%f\n",sum);
    V = (pow(2, m) / N) * sum - N;
    p_value = cephes_igamc((pow(2, m)-1) / 2, V / 2);

    fprintf(stats[TEST_POKERDETECT], "\t\t\t       POKERDETECT TEST\n");
	fprintf(stats[TEST_POKERDETECT], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_POKERDETECT], "\t\t COMPUTATIONAL INFORMATION:		  \n");
	fprintf(stats[TEST_POKERDETECT], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_POKERDETECT], "\t\t(a) Block length    (m) = %d\n", m);
	fprintf(stats[TEST_POKERDETECT], "\t\t(b) Sequence length (n) = %d\n", n);
	fprintf(stats[TEST_POKERDETECT], "\t\t(c) V                   = %f\n", V);
	fprintf(stats[TEST_POKERDETECT], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_POKERDETECT], "%s\t\tp_value = %f\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_POKERDETECT], "%f\n", p_value);
}