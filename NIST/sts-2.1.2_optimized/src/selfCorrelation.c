#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "../include/externs.h"

void
SelfCorrelation(int d, int n)
{
    // for(int i=0;i<n;i++){
    //     printf("%1d ",epsilon[i]);
    // }
    // printf("\n");

    int     i;
    double  V, n_d, sum, p_value, sqrt2 = 1.41421356237309504880;

    n_d = n - d;
    sum = 0.0;
    for (i = 0; i < n_d ; ++i) {
        sum += (epsilon[i] ^ epsilon[i + d]);
        // printf("%1d ",(epsilon[i] ^ epsilon[i + d]));
    }
    // printf("\n%f\n",sum);

    V = 2 * (sum - (n_d / 2)) / sqrt(n_d);
    p_value = erfc(fabs(V) / sqrt2);

	fprintf(stats[TEST_SELFCORRELATION], "\t\t\t      SELFCORRELATION TEST\n");
	fprintf(stats[TEST_SELFCORRELATION], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_SELFCORRELATION], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_SELFCORRELATION], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_SELFCORRELATION], "\t\t(a) shift length d      = %d\n", d);
	fprintf(stats[TEST_SELFCORRELATION], "\t\t(b) The nth partial sum = %d\n", (int)sum);
	fprintf(stats[TEST_SELFCORRELATION], "\t\t(c) V                   = %f\n", V);
	fprintf(stats[TEST_SELFCORRELATION], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_SELFCORRELATION], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_SELFCORRELATION]);
	fprintf(results[TEST_SELFCORRELATION], "%f\n", p_value); fflush(results[TEST_SELFCORRELATION]);
}

double sum1;
void
SelfCorrelation1(int d, int n)
{
    // for(int i=0;i<n;i++){
    //     printf("%1d ",epsilon[i]);
    // }
    // printf("\n");

    double  V, n_d, p_value, sqrt2 = 1.41421356237309504880;

    n_d = n - d;
    sum1 = 0.0;
    #pragma omp parallel for num_threads(8) reduction(+:sum1) 
    for (int i = 0; i < (int)n_d ; ++i) {
        sum1 += (epsilon[i] ^ epsilon[i + d]);
        // printf("%1d ",(epsilon[i] ^ epsilon[i + d]));
    }
    // printf("\n%f\n",sum1);

    V = 2 * (sum1 - (n_d / 2)) / sqrt(n_d);
    p_value = erfc(fabs(V) / sqrt2);

	fprintf(stats[TEST_SELFCORRELATION], "\t\t\t      SELFCORRELATION TEST\n");
	fprintf(stats[TEST_SELFCORRELATION], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_SELFCORRELATION], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_SELFCORRELATION], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_SELFCORRELATION], "\t\t(a) shift length d      = %d\n", d);
	fprintf(stats[TEST_SELFCORRELATION], "\t\t(b) The nth partial sum1 = %d\n", (int)sum1);
	fprintf(stats[TEST_SELFCORRELATION], "\t\t(c) V                   = %f\n", V);
	fprintf(stats[TEST_SELFCORRELATION], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_SELFCORRELATION], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_SELFCORRELATION]);
	fprintf(results[TEST_SELFCORRELATION], "%f\n", p_value); fflush(results[TEST_SELFCORRELATION]);
}