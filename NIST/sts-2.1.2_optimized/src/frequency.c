#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                          F R E Q U E N C Y  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
Frequency(int n)
{
	int		i;
	double	f, s_obs, p_value, sum, sqrt2 = 1.41421356237309504880;
	
	sum = 0.0;
	for ( i=0; i<n; i++ )
		sum += 2*(int)epsilon[i]-1;
	s_obs = fabs(sum)/sqrt(n);
	f = s_obs/sqrt2;
	p_value = erfc(f);

	fprintf(stats[TEST_FREQUENCY], "\t\t\t      FREQUENCY TEST\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t(a) The nth partial sum = %d\n", (int)sum);
	fprintf(stats[TEST_FREQUENCY], "\t\t(b) S_n/n               = %f\n", sum/n);
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_FREQUENCY]);
	fprintf(results[TEST_FREQUENCY], "%f\n", p_value); fflush(results[TEST_FREQUENCY]);
}

void
Frequency1(int n)
{
	int		i, j;
	double	f, s_obs, p_value, sum, sqrt2 = 1.41421356237309504880;
	
	sum = 0.0;
	int byteNum = n / 8; 
	int temp = 0;
	for (i = 0; i < byteNum; i++) {
		unsigned char local = byteEpsilon[i];
		temp += LU_byte_weight[local];
		// sum+=(LU_byte_weight[local] << 1) - 8;
	}
	sum = (temp << 1) - 8 * byteNum;
	if (n % 8 != 0) {
		unsigned char local = byteEpsilon[byteNum];
		for (j = 0; j < n - 8 * byteNum; j++) {
			sum += (((local >> (7-j)) & 1) << 1) - 1;
		}
	}
	
	s_obs = fabs(sum)/sqrt(n);
	f = s_obs/sqrt2;
	p_value = erfc(f);

	// write the result into file
	fprintf(stats[TEST_FREQUENCY], "\t\t\t      FREQUENCY TEST\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t(a) The nth partial sum = %d\n", (int)sum);
	fprintf(stats[TEST_FREQUENCY], "\t\t(b) S_n/n               = %f\n", sum/n);
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_FREQUENCY]);
	fprintf(results[TEST_FREQUENCY], "%f\n", p_value); fflush(results[TEST_FREQUENCY]);
}
