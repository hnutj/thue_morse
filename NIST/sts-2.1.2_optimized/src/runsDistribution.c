#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/cephes.h"
#include "../include/externs.h"

void
RunsDistribution(int n)
{
    // for(int i=0;i<n;i++){
    //     printf("%1d ",epsilon[i]);
    // }
    // printf("\n");
    int i, j, e, pos, k;
    unsigned char runFlag;
    double p_value, V, T;
    double *bi = NULL, *gi = NULL, *ei = NULL;

    // 1<=i<=n 时，随i增大，e单调递减
    for (i = 1; i <= n; ++i) {
        e = (double)(n - i + 3) / pow(2, i + 2);
        if (e >= 5) {
            k = i;
        }else{
            break;
        }
    }

	if ( ((bi = (double*)calloc(k+1, sizeof(double))) == NULL) || 
         ((gi = (double*)calloc(k+1, sizeof(double))) == NULL) ||
         ((ei = (double*)calloc(k+1, sizeof(double))) == NULL) ) {
		fprintf(stats[TEST_RUNSDISTRIBUTION], "\tInsufficient memory for required work space.\n");
		return;
	}

    runFlag = epsilon[0];
    j = 1;
    for (i = 1; i < n; ++i) {
        if (epsilon[i] != runFlag) {
            pos = j>k? k: j;
            if (runFlag == 0) {
                // printf("gi:pos=%d\n",pos);
                gi[pos] += 1;
            } else if (runFlag == 1) {
                // printf("bi:pos=%d\n",pos);
                bi[pos] += 1;
            }
            runFlag = epsilon[i];
            j = 1;
        } else {
            ++j;
        }
    }

    // record the last run
    pos = j>k? k: j;
    if (runFlag == 0) {
        // printf("gi:pos=%d\n",pos);
        gi[pos] += 1;
    } else if (runFlag == 1) {
        // printf("bi:pos=%d\n",pos);  
        bi[pos] += 1;
    }

    T = 0.0;
    for( i=1; i<=k; ++i ){
        T += bi[i]+gi[i];
    }

    for( i=1; i<k; ++i ){
        ei[i] = T / pow(2,i+1);
    }
    ei[k] = T / pow(2,k);

    V = 0.0;
    for (i = 1; i <= k; ++i) {
        double et = ei[i];
        V += pow(bi[i] - et, 2) / et;
        V += pow(gi[i] - et, 2) / et;
    }

    p_value = cephes_igamc(k - 1, V / 2);

	fprintf(stats[TEST_RUNSDISTRIBUTION], "\t\t\t      RUNSDISTRIBUTION TEST\n");
	fprintf(stats[TEST_RUNSDISTRIBUTION], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_RUNSDISTRIBUTION], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_RUNSDISTRIBUTION], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_RUNSDISTRIBUTION], "\t\t(a) k                   = %d\n", k);
	fprintf(stats[TEST_RUNSDISTRIBUTION], "\t\t(b) V                   = %f\n", V);
	fprintf(stats[TEST_RUNSDISTRIBUTION], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_RUNSDISTRIBUTION], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_RUNSDISTRIBUTION]);
	fprintf(results[TEST_RUNSDISTRIBUTION], "%f\n", p_value); fflush(results[TEST_RUNSDISTRIBUTION]);

    free(bi);
    free(gi);
    free(ei);
}