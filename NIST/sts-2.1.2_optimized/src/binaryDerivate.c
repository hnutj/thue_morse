#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <arpa/inet.h>
#include <omp.h>
#include "../include/externs.h"

void
BinaryDerivate(int k, int n)
{
    // for(int i=0;i<n;i++){
    //     printf("%1d ",epsilon[i]);
    // }
    // printf("\n");

    int i, j;
    double V, n_k, sum, p_value, sqrt2 = 1.41421356237309504880;
    BitSequence* sequence;

    sum = 0.0;
    n_k = n - k;

	if ( (sequence = (BitSequence *) malloc(n * sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}
    // 可前置到程序读入初始化部分
    memcpy(sequence, epsilon, n * sizeof(BitSequence));

    for (i = 0; i < k; ++i) {
        for (j = 0; j < n - i - 1; ++j) {
            sequence[j] = sequence[j] ^ sequence[j + 1];
        //     printf("%1d ",sequence[j]);           
        }
        // printf("\n");
    }
    for (i = 0; i < n_k; ++i) {
        sum += (2 * (int)sequence[i]) - 1;
    }
    // printf("\n%f\n",sum);

    V = sum / sqrt(n_k);
    p_value = erfc(fabs(V) / sqrt2);

	fprintf(stats[TEST_BINARYDERIVATE], "\t\t\t      BINARYDERIVATE TEST\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(a) power k             = %d\n", k);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(b) The nth partial sum = %d\n", (int)sum);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(c) V                   = %f\n", V);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_BINARYDERIVATE], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BINARYDERIVATE]);
	fprintf(results[TEST_BINARYDERIVATE], "%f\n", p_value); fflush(results[TEST_BINARYDERIVATE]);

	free(sequence);
}

void
BinaryDerivate1(int k, int n)
{
    // for(int i=0;i<n;i++){
    //     printf("%1d ",epsilon[i]);
    // }
    // printf("\n");
    
    int mask, m, endByte, endPos, oneNum;
    double V, n_k, sum, p_value, sqrt2 = 1.41421356237309504880;
    int * x;

    n_k = n - k;
	if ( (x = (int *) malloc((k+1) * sizeof(int))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}

	// 求杨辉三角的第k+1行，将系数向量存储在mask中
	x[0] = 1;
	mask = 1;
	for (int i = 1; i <= k; i++) {
		x[i] = x[i - 1] * (k - i + 1) / i;
		mask <<= 1;
		if (x[i] % 2 != 0) {
			mask |= 1;
		}
	}
	
    m = k + 1;
	endByte = (n - m) / 8;
	endPos = (n - m) % 8;
    oneNum = 0;	 // 对原始序列一次遍历,记录目标推导序列的‘1’的个数,以便直接求出序列和
	for (int i = 0; i < endByte; ++i) {
		for (int j = 0; j < 8; ++j) {
			int pos = 8 - j - m;
			int temp;
			if (pos >= 0) {
				// 涉及1byte
				temp = (byteEpsilon[i] >> pos) & mask;
			}
			else if (pos >= -8) {
				//涉及2bytes
				temp = ((byteEpsilon[i] << (-pos)) | (byteEpsilon[i + 1] >> (8 + pos))) & mask;
			}
			else if (pos >= -24) {
				//涉及3bytes
				temp = (htonl(((int*)(byteEpsilon + i))[0]) >> (24 + pos)) & mask;
			}
            // 求出推导序列中该比特值
			int count = 0;
			for (int p = 0; p < m; p++) {
				if (temp & 1) {
					count++;
				}
				temp >>= 1;
			}
			if (count & 1) {
				oneNum++;
				// printf("1 ");
			}
			// else {
			// 	printf("0 ");
			// }
		}
	}
	//处理endByte
	for (int j = 0; j <= endPos; ++j) {
		int pos = 8 - j - m;
		int temp;
		if (pos >= 0) {
			// 涉及1byte
			temp = (byteEpsilon[endByte] >> pos) & mask;
		}
		else if (pos >= -8) {
			//涉及2bytes
			temp = ((byteEpsilon[endByte] << (-pos)) | (byteEpsilon[endByte + 1] >> (8 + pos))) & mask;
		}
		else if (pos >= -24) {
			//涉及3bytes
			temp = (htonl(((int*)(byteEpsilon + endByte))[0]) >> (24 + pos)) & mask;
		}
        // 求出推导序列中该比特值
		int count = 0;
		for (int p = 0; p < m; p++) {
			if (temp & 1) {
				count++;
			}
			temp >>= 1;
		}
        if (count & 1) {
            oneNum++;
            // printf("1 ");
        }
        // else {
        //     printf("0 ");
        // }
	}
	sum = 2 * oneNum - n + k;

    V = sum / sqrt(n_k);
    p_value = erfc(fabs(V) / sqrt2);

	fprintf(stats[TEST_BINARYDERIVATE], "\t\t\t      BINARYDERIVATE TEST\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(a) power k             = %d\n", k);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(b) The nth partial sum = %d\n", (int)sum);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(c) V                   = %f\n", V);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_BINARYDERIVATE], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BINARYDERIVATE]);
	fprintf(results[TEST_BINARYDERIVATE], "%f\n", p_value); fflush(results[TEST_BINARYDERIVATE]);

    free(x);
}

int oneNum2;
void
BinaryDerivate2(int k, int n)
{
    // for(int i=0;i<n;i++){
    //     printf("%1d ",epsilon[i]);
    // }
    // printf("\n");
    
    int validNum;
    double V, n_k, sum, p_value, sqrt2 = 1.41421356237309504880;
    int * x, * mask;

    n_k = n - k;
	if ( ((x = (int *) malloc((k+1) * sizeof(int))) == NULL) ||
		 ((mask = (int *) malloc((k+1) * sizeof(int))) == NULL)) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}

	// 求杨辉三角的第k+1行，将有效位的系数向量存储在mask中
	x[0] = 1;
	mask[0] = 0;
	validNum = 1;
	for (int i = 1; i <= k; i++) {
		x[i] = x[i - 1] * (k - i + 1) / i;
		if (x[i] & 1) {
			mask[validNum++] = i;
		}
	}
	
    oneNum2 = 0;	 // 对原始序列一次遍历,记录目标推导序列的‘1’的个数,以便直接求出序列和
	#pragma omp parallel for reduction(+:oneNum2)
	for (int i = 0; i < (int)n_k; i++) {
		int count = 0;
		for (int j = 0; j < validNum; j++) {
			int pos = i + mask[j];
			if (epsilon[pos] & 1) {
			// if ((byteEpsilon[pos / 8] >> (7 - pos % 8)) & 1) {
				count++;
			}
		}
        if (count & 1) {
            oneNum2++;
            // printf("1 ");
        }
        // else {
        //     printf("0 ");
        // }
	}

	sum = 2 * oneNum2 - n + k;

    V = sum / sqrt(n_k);
    p_value = erfc(fabs(V) / sqrt2);

	fprintf(stats[TEST_BINARYDERIVATE], "\t\t\t      BINARYDERIVATE TEST\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(a) power k             = %d\n", k);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(b) The nth partial sum = %d\n", (int)sum);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(c) V                   = %f\n", V);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_BINARYDERIVATE], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BINARYDERIVATE]);
	fprintf(results[TEST_BINARYDERIVATE], "%f\n", p_value); fflush(results[TEST_BINARYDERIVATE]);

    free(x);
}

void
BinaryDerivate3(int k, int n)
{
    // for(int i=0;i<n;i++){
    //     printf("%1d ",epsilon[i]);
    // }
    // printf("\n");

    int i, j;
    double V, n_k, sum, p_value, sqrt2 = 1.41421356237309504880;
    BitSequence* sequence;

    sum = 0.0;
    n_k = n - k;

	if ( (sequence = (BitSequence *) malloc(n * sizeof(BitSequence))) == NULL ) {
			printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
			printf("Statistical Testing Aborted!\n");
			return;
		}

	// 若k为2的幂次减一，特殊处理提速
	if((k & (k+1)) == 0){
		sequence[0] = 0;
		for( i = 0; i <= k; ++i){
			sequence[0] ^= epsilon[i];
		}
		for (i = 1; i < n_k; ++i) {
			sequence[i] = sequence[i-1] ^ epsilon[i-1] ^ epsilon[i+k];
		}
	}
	else{
		memcpy(sequence, epsilon, n * sizeof(BitSequence));
		for (i = 0; i < k; ++i) {
			for (j = 0; j < n - i - 1; ++j) {
				sequence[j] = sequence[j] ^ sequence[j + 1];
			//     printf("%1d ",sequence[j]);           
			}
			// printf("\n");
		}
	}

	for (i = 0; i < n_k; ++i) {
		sum += (2 * (int)sequence[i]) - 1;
	}
	// printf("\n%f\n",sum);	
    V = sum / sqrt(n_k);
    p_value = erfc(fabs(V) / sqrt2);

	fprintf(stats[TEST_BINARYDERIVATE], "\t\t\t      BINARYDERIVATE TEST\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(a) power k             = %d\n", k);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(b) The nth partial sum = %d\n", (int)sum);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t(c) V                   = %f\n", V);
	fprintf(stats[TEST_BINARYDERIVATE], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_BINARYDERIVATE], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BINARYDERIVATE]);
	fprintf(results[TEST_BINARYDERIVATE], "%f\n", p_value); fflush(results[TEST_BINARYDERIVATE]);

	free(sequence);
}
