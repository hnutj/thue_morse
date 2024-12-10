#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <arpa/inet.h>
#include <omp.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  

extern void NonOverlappingHelper(int m, int n,int templateNum, unsigned char* byteEpsilon, unsigned short* templates, unsigned int*container);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          N O N O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
NonOverlappingTemplateMatchings(int m, int n)
{
	int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
						2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must 
	first be constructed, saved into files and then the corresponding 
	number of nonperiodic templates for that file be stored in the m-th 
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/
	unsigned int	bit, W_obs, nu[6], *Wj = NULL; 
	FILE			*fp = NULL;
	double			sum, chi2, p_value, lambda, pi[6], varWj;
	int				i, j, jj, k, match, SKIP, M, N, K = 5;
	char			directory[100];
	BitSequence		*sequence = NULL;

	N = 8;
	M = n/N;

	if ( (Wj = (unsigned int*)calloc(N, sizeof(unsigned int))) == NULL ) {
		fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
		fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
		return;
	}
	lambda = (M-m+1)/pow(2, m);
	varWj = M*(1.0/pow(2.0, m) - (2.0*m-1.0)/pow(2.0, 2.0*m));
	sprintf(directory, "templates/template%d", m);

	if ( ((isNegative(lambda)) || (isZero(lambda))) ||
		 ((fp = fopen(directory, "r")) == NULL) ||
		 ((sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == NULL) ) {
		fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
		fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
		fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
		fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
		if ( sequence != NULL )
			free(sequence);
	}
	else {
		fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
		fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");

		if ( numOfTemplates[m] < MAXNUMOFTEMPLATES )
			SKIP = 1;
		else
			SKIP = (int)(numOfTemplates[m]/MAXNUMOFTEMPLATES);
		numOfTemplates[m] = (int)numOfTemplates[m]/SKIP;
		
		sum = 0.0;
		for ( i=0; i<2; i++ ) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i];
		}
		pi[0] = sum;
		for ( i=2; i<=K; i++ ) {                      /* Compute Probabilities */
			pi[i-1] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i-1];
		}
		pi[K] = 1 - sum;

		for( jj=0; jj<MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]); jj++ ) {
			sum = 0;

			for ( k=0; k<m; k++ ) {
				fscanf(fp, "%d", &bit);
				sequence[k] = bit;
				fprintf(stats[TEST_NONPERIODIC], "%d", sequence[k]);
			}
			fprintf(stats[TEST_NONPERIODIC], " ");
			for ( k=0; k<=K; k++ )
				nu[k] = 0;
			for ( i=0; i<N; i++ ) {
				W_obs = 0;
				for ( j=0; j<M-m+1; j++ ) {
					match = 1;
					for ( k=0; k<m; k++ ) {
						if ( (int)sequence[k] != (int)epsilon[i*M+j+k] ) {
							match = 0;
							break;
						}
					}
					if ( match == 1 ) {
						W_obs++;
                        j += m-1;
                    }
				}
				Wj[i] = W_obs;
			}
			sum = 0;
			chi2 = 0.0;                                   /* Compute Chi Square */
			for ( i=0; i<N; i++ ) {
				if ( m == 10 )
					fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i]);
				else
					fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i]);
				chi2 += pow(((double)Wj[i] - lambda)/pow(varWj, 0.5), 2);
			}
			p_value = cephes_igamc(N/2.0, chi2/2.0);
		
			if ( isNegative(p_value) || isGreaterThanOne(p_value) )
				fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
			if ( SKIP > 1 )
				fseek(fp, (long)(SKIP-1)*2*m, SEEK_CUR);
			fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);
		}
	}
	
	fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);
	if ( sequence != NULL )
		free(sequence);

	free(Wj);
    if ( fp != NULL )
        fclose(fp);
}

void
NonOverlappingTemplateMatchings1(int m, int n)
{
	// m为待检测模板的长度,numOfTemplates记录这个模板长度下有多少个待检测的序列
	int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
						2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must 
	first be constructed, saved into files and then the corresponding 
	number of nonperiodic templates for that file be stored in the m-th 
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/

	unsigned int	bit, **Wj = NULL; 
	FILE			*fp = NULL;
	double			sum, chi2, p_value, lambda, pi[6], varWj;
	int				i, jj, SKIP, M, N, K = 5;
	char			directory[100];
	BitSequence		*sequence = NULL;

	// 锁定分块数为8
	N = 8;
	// 每个子块的长度
	M = n/N;

	// 均值
	lambda = (M-m+1)/pow(2, m);
	// 方差
	varWj = M*(1.0/pow(2.0, m) - (2.0*m-1.0)/pow(2.0, 2.0*m));
	// directory这个目录下存储的是这个m-bit长度的所有带检测模板
	sprintf(directory, "templates/template%d", m);

	if ( ((isNegative(lambda)) || (isZero(lambda))) ||
		 ((fp = fopen(directory, "r")) == NULL) ||
		 ((sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == NULL) ) {
		fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
		fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
		fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
		fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
		if ( sequence != NULL )
			free(sequence);
	}
	else {
		fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
		fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");

		if ( numOfTemplates[m] < MAXNUMOFTEMPLATES )
			SKIP = 1;
		else
			SKIP = (int)(numOfTemplates[m]/MAXNUMOFTEMPLATES);
		numOfTemplates[m] = (int)numOfTemplates[m]/SKIP;
		
		sum = 0.0;
		for ( i=0; i<2; i++ ) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i];
		}
		pi[0] = sum;
		for ( i=2; i<=K; i++ ) {                      /* Compute Probabilities */
			pi[i-1] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i-1];
		}
		pi[K] = 1 - sum;

		//预先读入
		int templateNum = MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]);
		unsigned short* templates;
		if ( (templates = (unsigned short*)calloc(templateNum, sizeof(unsigned short))) == NULL ) {
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
			fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
			return;
		}
		if ( (Wj = (unsigned int**)realloc(Wj, sizeof(unsigned int*)*templateNum)) == NULL ) {
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
			fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
			return;
		}
		for( jj=0; jj<templateNum; jj++){
			Wj[jj]=(unsigned int *)calloc(N,sizeof(unsigned int));
			templates[jj] = 0;
			for (int k = 0; k < m; k++) {
				fscanf(fp, "%1d", &bit);
				templates[jj] = (templates[jj] << 1) | bit;
			}
			if ( SKIP > 1 )
				fseek(fp, (long)(SKIP-1)*2*m, SEEK_CUR);
		}

		unsigned short mask = pow(2, m) - 1;

		// #pragma omp parallel for num_threads(8) private(j) 
		#pragma omp parallel for collapse(2)
		for( int jj=0; jj<templateNum; ++jj ) {
			for ( int j = 0; j < N; ++j) {
				unsigned short B = templates[jj]; //置于内层循环内外均可,置于内层为的是折叠循环便于omp并行
				//每个子块的索引区域：j*blockSize~(j+1)*blockSize-m
				int beginPos = j * M, endPos = beginPos + M - m;
				int beginByte = beginPos / 8, beginOff = beginPos % 8;
				int endByte = endPos / 8, endOff = endPos % 8;

				int W_obs = 0; //每个模板下每个子块有一个W_obs

				int p = beginByte, t=beginOff;
				while (p <= endByte) {
					int pos2 = (p == endByte ? endOff : 7);
					if( t > pos2){
						break;
					}
					while (t <= pos2) {
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
						if (test == B) {
							//match
							W_obs++;
							//int nextPos = p * 8 + t + m;
							p += (t + m) / 8;
							t = (t + m) % 8;
							break;
						}
						else if(t++ == pos2){
							t = 0;
							p++;
							break;
						}
					}
				}
				Wj[jj][j] = W_obs;
			}
		}

		for( jj=0; jj<templateNum; jj++ ) {
			fprintf(stats[TEST_NONPERIODIC], "%d", templates[jj]);
			fprintf(stats[TEST_NONPERIODIC], " ");
			chi2 = 0.0;                                   /* Compute Chi Square */
			for ( i=0; i<N; i++ ) {
				if ( m == 10 )
					fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[jj][i]);
				else
					fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[jj][i]);
				chi2 += pow(((double)Wj[jj][i] - lambda)/pow(varWj, 0.5), 2);
			}
			p_value = cephes_igamc(N/2.0, chi2/2.0);
		
			if ( isNegative(p_value) || isGreaterThanOne(p_value) )
				fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
			fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);
		}
	}

	fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);
	if ( sequence != NULL )
		free(sequence);

	free(Wj);
    if ( fp != NULL )
        fclose(fp);
}

void
NonOverlappingTemplateMatchings2(int m, int n)
{
	// m为待检测模板的长度,numOfTemplates记录这个模板长度下有多少个待检测的序列
	int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
						2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must 
	first be constructed, saved into files and then the corresponding 
	number of nonperiodic templates for that file be stored in the m-th 
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/

	unsigned int	bit; 
	FILE			*fp = NULL;
	double			sum, chi2, p_value, lambda, pi[6], varWj;
	int				i, jj, SKIP, M, N, K = 5;
	char			directory[100];
	BitSequence		*sequence = NULL;

	// 锁定分块数为8
	N = 8;
	// 每个子块的长度
	M = n/N;

	// 均值
	lambda = (M-m+1)/pow(2, m);
	// 方差
	varWj = M*(1.0/pow(2.0, m) - (2.0*m-1.0)/pow(2.0, 2.0*m));
	// directory这个目录下存储的是这个m-bit长度的所有带检测模板
	sprintf(directory, "templates/template%d", m);

	if ( ((isNegative(lambda)) || (isZero(lambda))) ||
		 ((fp = fopen(directory, "r")) == NULL) ||
		 ((sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == NULL) ) {
		fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
		fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
		fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
		fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
		if ( sequence != NULL )
			free(sequence);
	}
	else {
		fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
		fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");

		if ( numOfTemplates[m] < MAXNUMOFTEMPLATES )
			SKIP = 1;
		else
			SKIP = (int)(numOfTemplates[m]/MAXNUMOFTEMPLATES);
		numOfTemplates[m] = (int)numOfTemplates[m]/SKIP;
		
		sum = 0.0;
		for ( i=0; i<2; i++ ) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i];
		}
		pi[0] = sum;
		for ( i=2; i<=K; i++ ) {                      /* Compute Probabilities */
			pi[i-1] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i-1];
		}
		pi[K] = 1 - sum;

		//预先读入
		int templateNum = MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]);
		unsigned short* templates;
		if ( (templates = (unsigned short*)calloc(templateNum, sizeof(unsigned short))) == NULL ) {
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
			fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
			return;
		}
		for( jj=0; jj<templateNum; jj++){
			templates[jj] = 0;
			for (int k = 0; k < m; k++) {
				fscanf(fp, "%1d", &bit);
				templates[jj] = (templates[jj] << 1) | bit;
			}
			if ( SKIP > 1 )
				fseek(fp, (long)(SKIP-1)*2*m, SEEK_CUR);
		}

		unsigned int  *container = NULL;
		container = (unsigned int*)calloc(templateNum * N, sizeof(unsigned int));

		NonOverlappingHelper(m, n, templateNum, byteEpsilon, templates, container);

		for( jj=0; jj<templateNum; jj++ ) {
			fprintf(stats[TEST_NONPERIODIC], "%d", templates[jj]);
			fprintf(stats[TEST_NONPERIODIC], " ");
			chi2 = 0.0;                                   /* Compute Chi Square */
			for ( i=0; i<N; i++ ) {
				if ( m == 10 )
					fprintf(stats[TEST_NONPERIODIC], "%3d  ", container[jj*N+i]);
				else
					fprintf(stats[TEST_NONPERIODIC], "%4d ", container[jj*N+i]);
				chi2 += pow(((double)container[jj*N+i] - lambda)/pow(varWj, 0.5), 2);
			}
			p_value = cephes_igamc(N/2.0, chi2/2.0);
		
			if ( isNegative(p_value) || isGreaterThanOne(p_value) )
				fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
			fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);
		}
		free(container);
	}

	fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);
	if ( sequence != NULL )
		free(sequence);

	if ( fp != NULL )
        fclose(fp);
}