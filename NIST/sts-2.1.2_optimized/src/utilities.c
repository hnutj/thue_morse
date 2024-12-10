/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
U T I L I T I E S
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/generators.h"
#include "../include/stat_fncs.h"

int testJumpHop = 0;
int testDM = 0;

int
displayGeneratorOptions()
{
	int		option = 0;

	printf("           G E N E R A T O R    S E L E C T I O N \n");
	printf("           ______________________________________\n\n");
	printf("    [0] Input File                 [1] Linear Congruential\n");
	printf("    [2] Quadratic Congruential I   [3] Quadratic Congruential II\n");
	printf("    [4] Cubic Congruential         [5] XOR\n");
	printf("    [6] Modular Exponentiation     [7] Blum-Blum-Shub\n");
	printf("    [8] Micali-Schnorr             [9] G Using SHA-1\n\n");
	printf("   Enter Choice: ");
	scanf("%d", &option);
	printf("\n\n");

	return option;
}

// Input the name of the PRNG (0 indicates the user input file)
int
generatorOptions(char** streamFile)
{
	char	file[200];
	int		option = NUMOFGENERATORS+1;
	FILE	*fp;
	
	while ( (option < 0) || (option > NUMOFGENERATORS) ) {
		option = displayGeneratorOptions();
		switch( option ) {
			case 0:
				printf("\t\tUser Prescribed Input File: ");
				scanf("%s", file);
				*streamFile = (char*)calloc(200, sizeof(char));
				sprintf(*streamFile, "%s", file);
				printf("\n");
				if ( (fp = fopen(*streamFile, "r")) == NULL ) {
					printf("File Error:  file %s could not be opened.\n",  *streamFile);
					exit(-1);
				}
				else
					fclose(fp);
				break;
			case 1:
				*streamFile = "Linear-Congruential";
				break;
			case 2:
				*streamFile = "Quadratic-Congruential-1";
				break;
			case 3:
				*streamFile = "Quadratic-Congruential-2";
				break;
			case 4:
				*streamFile = "Cubic-Congruential";
				break;
			case 5:
				*streamFile = "XOR";
				break;
			case 6:
				*streamFile = "Modular-Exponentiation";
				break;
			case 7:
				*streamFile = "Blum-Blum-Shub";
				break;
			case 8:
				*streamFile = "Micali-Schnorr";
				break;
			case 9:
				*streamFile = "G using SHA-1";
				break;
				
			/* INTRODUCE NEW PRNG NAMES HERE */
			/*
			case 10:  *streamFile = "myNewPRNG";
				break;
			*/
			default:
				printf("Error:  Out of range - Try again!\n");
				break;
		}
	}
	return option;
}

// Select the test (use testVector to store, 1 means execute the test)
void
chooseTests()
{
	int		i;
	
	printf("                S T A T I S T I C A L   T E S T S\n");
	printf("                _________________________________\n\n");
	printf("    [01] Frequency                       [02] Block Frequency\n");
	printf("    [03] Cumulative Sums                 [04] Runs\n");
	printf("    [05] Longest Run of Ones             [06] Rank\n");
	printf("    [07] Discrete Fourier Transform      [08] Nonperiodic Template Matchings\n");
	printf("    [09] Overlapping Template Matchings  [10] Universal Statistical\n");
	printf("    [11] Approximate Entropy             [12] Random Excursions\n");
	printf("    [13] Random Excursions Variant       [14] Serial\n");
	printf("    [15] Linear Complexity               [16] Deviation Measure I\n");
    printf("    [17] Deviation Measure II            [18] Deviation Measure I CS\n");
    printf("    [19] Deviation Measure II CS         [20] Odd Hop Sum CS\n");
    printf("    [21] Even Hop Sum CS                 [22] Jump Complexity CS\n");
    printf("    [23] SelfCorrelation                 [24] BinaryDerivate\n");
    printf("    [25] PokerDetect                     [26] RunsDistribution\n");
	printf("         INSTRUCTIONS\n");
	printf("            Enter 0 if you DO NOT want to apply all of the\n");
	printf("            statistical tests to each sequence and 1 if you DO.\n\n");
	printf("   Enter Choice: ");
	scanf("%d", &testVector[0]);
	printf("\n");
	if ( testVector[0] == 1 )
		for( i=1; i<=NUMOFTESTS; i++ )
			testVector[i] = 1;
	else {
		printf("         INSTRUCTIONS\n");
		printf("            Enter a 0 or 1 to indicate whether or not the numbered statistical\n");
		printf("            test should be applied to each sequence.\n\n");
		printf("      12345678911111111112222222\n");
		printf("               01234567890123456\n");
		printf("      ");
		for ( i=1; i<=NUMOFTESTS; i++ ) 
			scanf("%1d", &testVector[i]);
		printf("\n\n");
	}
}

// adjust parameters
void
fixParameters()
{
	int		counter, testid;
	
	//  Check to see if any parameterized tests are selected
	if ( (testVector[TEST_BLOCK_FREQUENCY] != 1) && (testVector[TEST_NONPERIODIC] != 1) && 
		 (testVector[TEST_OVERLAPPING] != 1) && (testVector[TEST_APEN] != 1) &&
		 (testVector[TEST_SERIAL] != 1) && (testVector[TEST_LINEARCOMPLEXITY] != 1) &&
		 (testVector[TEST_SELFCORRELATION] != 1) && (testVector[TEST_BINARYDERIVATE] != 1) &&
		 (testVector[TEST_POKERDETECT] != 1) )
			return;
		
	do {
		counter = 1;
		printf("        P a r a m e t e r   A d j u s t m e n t s\n");
		printf("        -----------------------------------------\n");
		if ( testVector[TEST_BLOCK_FREQUENCY] == 1 )
			printf("    [%d] Block Frequency Test - block length(M):         %d\n", counter++, tp.blockFrequencyBlockLength);
		if ( testVector[TEST_NONPERIODIC] == 1 )
			printf("    [%d] NonOverlapping Template Test - block length(m): %d\n", counter++, tp.nonOverlappingTemplateBlockLength);
		if ( testVector[TEST_OVERLAPPING] == 1 )
			printf("    [%d] Overlapping Template Test - block length(m):    %d\n", counter++, tp.overlappingTemplateBlockLength);
		if ( testVector[TEST_APEN] == 1 )
			printf("    [%d] Approximate Entropy Test - block length(m):     %d\n", counter++, tp.approximateEntropyBlockLength);
		if ( testVector[TEST_SERIAL] == 1 )
			printf("    [%d] Serial Test - block length(m):                  %d\n", counter++, tp.serialBlockLength);
		if ( testVector[TEST_LINEARCOMPLEXITY] == 1 )
			printf("    [%d] Linear Complexity Test - block length(M):       %d\n", counter++, tp.linearComplexitySequenceLength);
		if ( testVector[TEST_SELFCORRELATION] == 1 )
			printf("    [%d] SelfCorrelation Test - shift length(d):         %d\n", counter++, tp.selfCorrelationShiftLength);
		if ( testVector[TEST_BINARYDERIVATE] == 1 )
			printf("    [%d] binaryDerivate Test - power(k):                 %d\n", counter++, tp.binaryDerivatePower);
		if ( testVector[TEST_POKERDETECT] == 1 )
			printf("    [%d] PokerDetect Test - block length(m):             %d\n", counter++, tp.pokerDetectBlockLength);
		printf("\n");
		printf("   Select Test (0 to continue): ");
		scanf("%1d", &testid);
		printf("\n");
		
		counter = 0;
		if ( testVector[TEST_BLOCK_FREQUENCY] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Block Frequency Test block length: ");
				scanf("%d", &tp.blockFrequencyBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_NONPERIODIC] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter NonOverlapping Template Test block Length: ");
				scanf("%d", &tp.nonOverlappingTemplateBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_OVERLAPPING] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Overlapping Template Test block Length: ");
				scanf("%d", &tp.overlappingTemplateBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_APEN] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Approximate Entropy Test block Length: ");
				scanf("%d", &tp.approximateEntropyBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_SERIAL] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Serial Test block Length: ");
				scanf("%d", &tp.serialBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_LINEARCOMPLEXITY] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Linear Complexity Test block Length: ");
				scanf("%d", &tp.linearComplexitySequenceLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_SELFCORRELATION] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter SelfCorrelation Test Shift Length: ");
				scanf("%d", &tp.selfCorrelationShiftLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_BINARYDERIVATE] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter binaryDerivate Test Power: ");
				scanf("%d", &tp.binaryDerivatePower);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_POKERDETECT] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter PokerDetect Block Length: ");
				scanf("%d", &tp.pokerDetectBlockLength);
				printf("\n");
				continue;
			}
		}

	} while ( testid != 0 );
}

// Input File Format if user input the PRNGs
void
fileBasedBitStreams(char *streamFile)
{
	FILE	*fp;
	int		mode;
	
	printf("   Input File Format:\n");
	printf("    [0] ASCII - A sequence of ASCII 0's and 1's\n");
	printf("    [1] Binary - Each byte in data file contains 8 bits of data\n\n");
	printf("   Select input mode:  ");
	scanf("%1d", &mode);
	printf("\n");
	if ( mode == 0 ) {
		// open the file as read-only
		if ( (fp = fopen(streamFile, "r")) == NULL ) {
			printf("ERROR IN FUNCTION fileBasedBitStreams:  file %s could not be opened.\n",  streamFile);
			exit(-1);
		}
		readBinaryDigitsInByteAndASCIIFormat(fp, streamFile);
		// readBinaryDigitsInASCIIFormat(fp, streamFile);
		fclose(fp);
	}
	else if ( mode == 1 ) {
		// open the binary file as read-only
		if ( (fp = fopen(streamFile, "rb")) == NULL ) {
			printf("ERROR IN FUNCTION fileBasedBitStreams:  file %s could not be opened.\n", streamFile);
			exit(-1);
		}
		readHexDigitsInByteAndBinaryFormat(fp);
		// readHexDigitsInBinaryFormat(fp);
		fclose(fp);
	}
}

// read the file with a sequence of ASCII 0's and 1's
void
readBinaryDigitsInASCIIFormat(FILE *fp, char *streamFile)
{
	int		i, j, num_0s, num_1s, bitsRead, bit;
	
	// alloc the space to epsilon
	if ( (epsilon = (BitSequence *) calloc(tp.n, sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}
	printf("     Statistical Testing In Progress.........\n\n");   
	for ( i=0; i<tp.numOfBitStreams; i++ ) {
		num_0s = 0;
		num_1s = 0;
		bitsRead = 0;
		for ( j=0; j<tp.n; j++ ) {
			if ( fscanf(fp, "%1d", &bit) == EOF ) {
				printf("ERROR:  Insufficient data in file %s.  %d bits were read.\n", streamFile, bitsRead);
				fclose(fp);
				free(epsilon);
				return;
			}
			else {
				bitsRead++;
				if ( bit == 0 ) 
					num_0s++;
				else 
					num_1s++;
				epsilon[j] = bit;
			}
		}
		// write the data into freqfp.txt
		fprintf(freqfp, "\t\tBITSREAD = %d 0s = %d 1s = %d\n", bitsRead, num_0s, num_1s);
		nist_test_suite();
	}
	free(epsilon);
}

void
readBinaryDigitsInByteAndASCIIFormat(FILE *fp, char *streamFile)
{
	for(int i=0;i<NUMOFTESTS;i++){
		timeCounter[i] = 0.0;
	}
	int		i, j, num_0s, num_1s, bitsRead, bit;
	int     byteNum = tp.n / 8; //能这样舍去尾巴吗？
	int     byteTotalNum = tp.n % 8 == 0 ? byteNum + 4 : byteNum + 5;

	// alloc the space to epsilon
	if ( (epsilon = (BitSequence *) calloc(tp.n, sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}
	if ( (byteEpsilon = (BitSequence *) calloc(byteTotalNum, sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}

	printf("     Statistical Testing In Progress.........\n\n");   
	for ( i=0; i<tp.numOfBitStreams; i++ ) {
		num_0s = 0;
		num_1s = 0;
		bitsRead = 0;
		int byteCount = -1;
		long off1 = 0, off2 = 0;

		// 记录此时文件偏移
		off1 = ftell(fp);

		for ( j=0; j<tp.n; j++ ) {
			// 计算当前字节编号（从0开始）
			if( j%8 == 0 && byteCount < byteNum){
				byteCount++;
				byteEpsilon[byteCount] = 0;
			}
			
			if ( fscanf(fp, "%1d", &bit) == EOF ) {
				printf("ERROR:  Insufficient data in file %s.  %d bits were read.\n", streamFile, bitsRead);
				fclose(fp);
				free(epsilon);
				return;
			}
			else {
				bitsRead++;
				if ( bit == 0 ) 
					num_0s++;
				else 
					num_1s++;
				epsilon[j] = bit;
				byteEpsilon[byteCount] = (byteEpsilon[byteCount] << 1) | bit;
			}
		}
		// 记录此时文件偏移
		off2 = ftell(fp);
		byteEpsilon[byteCount+1]=byteEpsilon[byteCount+2]=byteEpsilon[byteCount+3]=byteEpsilon[byteCount+4]=0;

		int begin = byteNum;
		int left = tp.n - 8 * byteNum;
		// 跳回开头
		fseek(fp, off1, SEEK_SET);	
		for ( int k = 0; k < 31; k++) {
			if (fscanf(fp, "%1d", &bit) == EOF) {
				printf("ERROR:  Insufficient data in file %s.\n", streamFile);
				fclose(fp);
				return ;
			}else {
				byteEpsilon[begin] = (byteEpsilon[begin] << 1) | bit;
			}
			if (++left % 8 == 0) {
				begin++;
			}
		}

		//补齐移位
		int left1 = 8 * byteTotalNum - (tp.n + 31);
		byteEpsilon[byteTotalNum-1] = byteEpsilon[byteTotalNum - 1] << left1;

		// write the data into freqfp.txt
		fprintf(freqfp, "\t\tBITSREAD = %d 0s = %d 1s = %d\n", bitsRead, num_0s, num_1s);
		printf("bitStream:%d\n",i);
		// 恢复文件偏移
		fseek(fp, off2, SEEK_SET);
		nist_test_suite();
		if( (i+1)%20 == 0 || i==tp.numOfBitStreams-1){
			printf("the %dth bitStream has completed the test\n",i+1);
		}
	}
	printTime();
	free(epsilon);
}

void
readHexDigitsInByteAndBinaryFormat(FILE *fp)
{
	for(int i=0;i<NUMOFTESTS;i++){
		timeCounter[i] = 0.0;
	}

	int		i, j, k, done, bit, num_0s, num_1s, bitsRead;
	BYTE	buffer[4];	
	int     byteNum = tp.n / 8; 
	int     byteTotalNum = tp.n % 8 == 0 ? byteNum + 4 : byteNum + 5;

	// alloc the space to epsilon
	if ( (epsilon = (BitSequence *) calloc(tp.n, sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}
	if ( (byteEpsilon = (BitSequence *) calloc(byteTotalNum, sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}

	printf("     Statistical Testing In Progress.........\n\n");   
	for ( i=0; i<tp.numOfBitStreams; i++ ) {
		num_0s = 0;
		num_1s = 0;
		bitsRead = 0;
		done = 0;
		int byteCount = 0;
		int threshold = 0;
		long off1 = 0, off2 = 0;

		// 记录此时文件偏移
		off1 = ftell(fp);

		do {
			if ( fread(buffer, sizeof(unsigned char), 4, fp) != 4 ) {
				printf("READ ERROR:  Insufficient data in file.\n");
				free(epsilon);
				return;
			}
			// done = convertToBits(buffer, 32, tp.n, &num_0s, &num_1s, &bitsRead);
			for ( k=0; k<4 ; k++ ) {
				int mask = 0x80;
				if(byteCount == byteNum){
					threshold = 1;
					byteEpsilon[byteCount] = 0;
				}else{
					byteEpsilon[byteCount++] = *(buffer + k) ;
				}

				for ( j=0; j<8; j++ ) {
					if ( *(buffer + k) & mask ) {
						bit = 1;
						num_1s++;
					}
					else {
						bit = 0;
						num_0s++;
					}
					mask >>= 1;
					epsilon[bitsRead++] = bit;
					if(threshold){
						byteEpsilon[byteCount] = (byteEpsilon[byteCount] << 1) | bit;
					}
					if ( bitsRead == tp.n ){
						done = 1;
						break;	
					}
				}
				if(done){
					break;
				}
			}
		} while ( !done );
		fprintf(freqfp, "\t\tBITSREAD = %d 0s = %d 1s = %d\n", bitsRead, num_0s, num_1s);
		
		// 记录此时文件偏移
		off2 = ftell(fp);
		if( tp.n % 8 == 0) byteCount--;
		byteEpsilon[byteCount+1]=byteEpsilon[byteCount+2]=byteEpsilon[byteCount+3]=byteEpsilon[byteCount+4]=0;
		int begin = byteNum;
		int left = tp.n - 8 * byteNum;

		// 跳回开头
		fseek(fp, off1, SEEK_SET);	
		fread(buffer, sizeof(unsigned char), 4, fp);
		for ( k=0; k<4 ; k++ ) {
			int mask = 0x80;
			for ( j=0; j<8; j++ ) {
				if(k == 3 && j == 7) break;
				if ( *(buffer + k) & mask ) {
					bit = 1;
				}
				else {
					bit = 0;
				}
				mask >>= 1;
				byteEpsilon[begin] = (byteEpsilon[begin] << 1) | bit;
				if (++left % 8 == 0) {
					begin++;
				}
			}
		}

		//补齐移位
		int left1 = 8 * byteTotalNum - (tp.n + 31);
		byteEpsilon[byteTotalNum-1] = byteEpsilon[byteTotalNum - 1] << left1;

		// 恢复文件偏移
		fseek(fp, off2, SEEK_SET);
		nist_test_suite();
		if( (i+1)%20 == 0 || i==tp.numOfBitStreams-1){
			printf("the %dth bitStream has completed the test\n",i+1);
		}
	}
	printTime();
	free(epsilon);
}

// read the file and each byte in the file contains 8 bits of data
void
readHexDigitsInBinaryFormat(FILE *fp)
{
	int		i, done, num_0s, num_1s, bitsRead;
	BYTE	buffer[4];
	
	if ( (epsilon = (BitSequence *) calloc(tp.n,sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		return;
	}

	printf("     Statistical Testing In Progress.........\n\n");   
	for ( i=0; i<tp.numOfBitStreams; i++ ) {
		num_0s = 0;
		num_1s = 0;
		bitsRead = 0;
		done = 0;
		do {
			if ( fread(buffer, sizeof(unsigned char), 4, fp) != 4 ) {
				printf("READ ERROR:  Insufficient data in file.\n");
				free(epsilon);
				return;
			}
			done = convertToBits(buffer, 32, tp.n, &num_0s, &num_1s, &bitsRead);
		} while ( !done );
		fprintf(freqfp, "\t\tBITSREAD = %d 0s = %d 1s = %d\n", bitsRead, num_0s, num_1s);
		
		nist_test_suite();
		
	}
	free(epsilon);
}


int
convertToBits(BYTE *x, int xBitLength, int bitsNeeded, int *num_0s, int *num_1s, int *bitsRead)
{
	int		i, j, count, bit;
	BYTE	mask;
	int		zeros, ones;

	count = 0;
	zeros = ones = 0;
	for ( i=0; i<(xBitLength+7)/8; i++ ) {
		mask = 0x80;
		for ( j=0; j<8; j++ ) {
			if ( *(x+i) & mask ) {
				bit = 1;
				(*num_1s)++;
				ones++;
			}
			else {
				bit = 0;
				(*num_0s)++;
				zeros++;
			}
			mask >>= 1;
			epsilon[*bitsRead] = bit;
			(*bitsRead)++;
			if ( *bitsRead == bitsNeeded )
				return 1;
			if ( ++count == xBitLength )
				return 0;
		}
	}
	
	return 0;
}

// open the output file and input the number of bitstreams
void
openOutputStreams(int option)
{
	int		i, numOfBitStreams, numOfOpenFiles = 0;
	char	freqfn[200], summaryfn[200], statsDir[200], resultsDir[200];
	
	// There must be freq.txt and finalAnalysisReport.txt
	sprintf(freqfn, "experiments/%s/freq.txt", generatorDir[option]);
	if ( (freqfp = fopen(freqfn, "w")) == NULL ) {
		printf("\t\tMAIN:  Could not open freq file: <%s>", freqfn);
		exit(-1);
	}
	sprintf(summaryfn, "experiments/%s/finalAnalysisReport.txt", generatorDir[option]);
	if ( (summary = fopen(summaryfn, "w")) == NULL ) {
		printf("\t\tMAIN:  Could not open stats file: <%s>", summaryfn);
		exit(-1);
	}
	
	// There will be some txt in the chosen test dir 
	for( i=1; i<=NUMOFTESTS; i++ ) {
		if ( testVector[i] == 1 ) {
			sprintf(statsDir, "experiments/%s/%s/stats.txt", generatorDir[option], testNames[i]);
			sprintf(resultsDir, "experiments/%s/%s/results.txt", generatorDir[option], testNames[i]);
			if ( (stats[i] = fopen(statsDir, "w")) == NULL ) {	/* STATISTICS LOG */
				printf("ERROR: LOG FILES COULD NOT BE OPENED.\n");
				printf("       MAX # OF OPENED FILES HAS BEEN REACHED = %d\n", numOfOpenFiles);
				printf("-OR-   THE OUTPUT DIRECTORY DOES NOT EXIST.\n");
				exit(-1);
			}
			else
				numOfOpenFiles++;
			if ( (results[i] = fopen(resultsDir, "w")) == NULL ) {	/* P_VALUES LOG   */
				 printf("ERROR: LOG FILES COULD NOT BE OPENED.\n");
				 printf("       MAX # OF OPENED FILES HAS BEEN REACHED = %d\n", numOfOpenFiles);
				 printf("-OR-   THE OUTPUT DIRECTORY DOES NOT EXIST.\n");
				 exit(-1);
			}
			else
				numOfOpenFiles++;
		}
	}
	printf("   How many bitstreams? ");
	scanf("%d", &numOfBitStreams);
	tp.numOfBitStreams = numOfBitStreams;
	printf("\n");
}

// Generate(input) random sequences
void
invokeTestSuite(int option, char *streamFile)
{
	fprintf(freqfp, "________________________________________________________________________________\n\n");
	fprintf(freqfp, "\t\tFILE = %s\t\tALPHA = %6.4f\n", streamFile, ALPHA);
	fprintf(freqfp, "________________________________________________________________________________\n\n");
	if ( option != 0 )
		printf("     Statistical Testing In Progress.........\n\n");
	switch( option ) {
		case 0:
			fileBasedBitStreams(streamFile);
			break;
		case 1:
			lcg();
			break;
		case 2:
			quadRes1();
			break;
		case 3:
			quadRes2();
			break;
		case 4:
			cubicRes();
			break;
		case 5:
			exclusiveOR();
			break;
		case 6:
			modExp();
			break;
		case 7:
			bbs();
			break;
		case 8:
			micali_schnorr();
			break;
		case 9:
			SHA1();
			break;
			
		/* INTRODUCE NEW PSEUDO RANDOM NUMBER GENERATORS HERE */
			
		default:
			printf("Error in invokeTestSuite!\n");
			break;
	}
	printf("     Statistical Testing Complete!!!!!!!!!!!!\n\n");
}

void
nist_test_suite()
{
	struct timeval head, tail; 
	if ( (testVector[0] == 1) || (testVector[TEST_FREQUENCY] == 1) ) {
		gettimeofday(&head, NULL);		
		// Frequency(tp.n);
		Frequency1(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[0] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_BLOCK_FREQUENCY] == 1) ){
		gettimeofday(&head, NULL);		
		// BlockFrequency(tp.blockFrequencyBlockLength, tp.n);
		BlockFrequency1(tp.blockFrequencyBlockLength, tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[1] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	} 
	
	if ( (testVector[0] == 1) || (testVector[TEST_CUSUM] == 1) ){
		gettimeofday(&head, NULL);		
		// CumulativeSums(tp.n);
		CumulativeSums1(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[2] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_RUNS] == 1) ){
		gettimeofday(&head, NULL);
		// Runs(tp.n); 
		Runs1(tp.n); 
		gettimeofday(&tail, NULL);
		timeCounter[3] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_LONGEST_RUN] == 1) ){
		gettimeofday(&head, NULL);	
		// LongestRunOfOnes(tp.n);
		LongestRunOfOnes1(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[4] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_RANK] == 1) ){
		gettimeofday(&head, NULL);		
		// Rank(tp.n);
		Rank1(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[5] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_FFT] == 1) ){
		gettimeofday(&head, NULL);		
		// DiscreteFourierTransform(tp.n);
		// DiscreteFourierTransform1(tp.n);
		DiscreteFourierTransform2(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[6] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_NONPERIODIC] == 1) ){
		gettimeofday(&head, NULL);		
		// NonOverlappingTemplateMatchings(tp.nonOverlappingTemplateBlockLength, tp.n);
		NonOverlappingTemplateMatchings2(tp.nonOverlappingTemplateBlockLength, tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[7] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_OVERLAPPING] == 1) ){
		gettimeofday(&head, NULL);		
		// OverlappingTemplateMatchings(tp.overlappingTemplateBlockLength, tp.n);
		OverlappingTemplateMatchings1(tp.overlappingTemplateBlockLength, tp.n);		
		gettimeofday(&tail, NULL);
		timeCounter[8] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_UNIVERSAL] == 1) ){
		gettimeofday(&head, NULL);		
		// Universal(tp.n);
		Universal1(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[9] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_APEN] == 1) ){
		gettimeofday(&head, NULL);
		// ApproximateEntropy(tp.approximateEntropyBlockLength, tp.n);
		ApproximateEntropy1(tp.approximateEntropyBlockLength, tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[10] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_RND_EXCURSION] == 1) ){
		gettimeofday(&head, NULL);		
		// RandomExcursions(tp.n);
		RandomExcursions1(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[11] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_RND_EXCURSION_VAR] == 1) ){
		gettimeofday(&head, NULL);
		// RandomExcursionsVariant(tp.n);
		RandomExcursionsVariant1(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[12] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_SERIAL] == 1) ){
		gettimeofday(&head, NULL);
		// Serial(tp.serialBlockLength,tp.n);
		Serial1(tp.serialBlockLength,tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[13] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	
	if ( (testVector[0] == 1) || (testVector[TEST_LINEARCOMPLEXITY] == 1) ){
		gettimeofday(&head, NULL);
		// LinearComplexity(tp.linearComplexitySequenceLength, tp.n);
		LinearComplexity1(tp.linearComplexitySequenceLength, tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[14] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}

	if ((testVector[0] == 1) || 
	((testVector[TEST_DEVIATIONMEARSURE1] == 1) && (testVector[TEST_DEVIATIONMEARSURE2] == 1) &&
	 (testVector[TEST_DEVIATIONMEARSURE1CS] == 1) && (testVector[TEST_DEVIATIONMEARSURE2CS] == 1))){
		testDM = 1;
		gettimeofday(&head, NULL);
        DeviationMeasureTotal1(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[15] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	else{
		if ( testVector[TEST_DEVIATIONMEARSURE1] == 1 ){
			gettimeofday(&head, NULL);
			// DeviationMeasure1(tp.n);
			DeviationMeasure1_1(tp.n);
			gettimeofday(&tail, NULL);
			timeCounter[15] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
		}
			
		if ( testVector[TEST_DEVIATIONMEARSURE2] == 1 ){
			gettimeofday(&head, NULL);
			// DeviationMeasure2(tp.n);
			DeviationMeasure2_1(tp.n);
			gettimeofday(&tail, NULL);
			timeCounter[16] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
		}	
		
		if ( testVector[TEST_DEVIATIONMEARSURE1CS] == 1 ){
			gettimeofday(&head, NULL);
			// DeviationMeasure1CS(tp.n);
			DeviationMeasure1CS1(tp.n);
			gettimeofday(&tail, NULL);
			timeCounter[17] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
		}	
		
		if ( testVector[TEST_DEVIATIONMEARSURE2CS] == 1 ){
			gettimeofday(&head, NULL);
			// DeviationMeasure2CS(tp.n);
			DeviationMeasure2CS1(tp.n);
			gettimeofday(&tail, NULL);
			timeCounter[18] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
		}		
	}

	if ((testVector[0] == 1) || 
	((testVector[TEST_ODDHOPSUM] == 1) && (testVector[TEST_EVENHOPSUM] == 1) && (testVector[TEST_JUMP] == 1))){
		testJumpHop = 1;
		gettimeofday(&head, NULL);
        JumpHopComplexity1(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[19] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
	else{
		if (testVector[TEST_ODDHOPSUM] == 1){
			gettimeofday(&head, NULL);
			// Oddhopsum(tp.n);
			Oddhopsum1(tp.n);
			gettimeofday(&tail, NULL);
			timeCounter[19] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
		}

		if (testVector[TEST_EVENHOPSUM] == 1){
			gettimeofday(&head, NULL);
			// Evenhopsum(tp.n);
			Evenhopsum1(tp.n);
			gettimeofday(&tail, NULL);
			timeCounter[20] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
		}

		if (testVector[TEST_JUMP] == 1){
			gettimeofday(&head, NULL);
			// Jump(tp.n);
			Jump1(tp.n);
			gettimeofday(&tail, NULL);
			timeCounter[21] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
		}
	}

	if ( (testVector[0] == 1) || (testVector[TEST_SELFCORRELATION] == 1) ){
		gettimeofday(&head, NULL);
		// SelfCorrelation(tp.selfCorrelationShiftLength, tp.n);
		SelfCorrelation1(tp.selfCorrelationShiftLength, tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[22] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}

	if ( (testVector[0] == 1) || (testVector[TEST_BINARYDERIVATE] == 1) ){
		gettimeofday(&head, NULL);
		// BinaryDerivate(tp.binaryDerivatePower, tp.n);
		BinaryDerivate2(tp.binaryDerivatePower, tp.n);		
		gettimeofday(&tail, NULL);
		timeCounter[23] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}

	if ( (testVector[0] == 1) || (testVector[TEST_POKERDETECT] == 1) ){
		gettimeofday(&head, NULL);
		// PokerDetect(tp.pokerDetectBlockLength, tp.n);
		PokerDetect1(tp.pokerDetectBlockLength, tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[24] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}

	if ( (testVector[0] == 1) || (testVector[TEST_RUNSDISTRIBUTION] == 1) ){
		gettimeofday(&head, NULL);
		RunsDistribution(tp.n);
		// RunsDistribution1(tp.n);
		gettimeofday(&tail, NULL);
		timeCounter[25] += (tail.tv_sec-head.tv_sec)*1000.0+(tail.tv_usec-head.tv_usec)/1000.0;
	}
}

void
printTime()
{
	//化为秒单位
	double sum = 0;
	for(int i=0;i<NUMOFTESTS;i++){
		timeCounter[i]/=1000.0;
		sum += timeCounter[i];
	}

	printf("\nFrequency Test: time(s) %lf\n",timeCounter[0]);
	printf("BlockFrequency Test: time(s) %lf\n",timeCounter[1]);
	printf("CumulativeSums Test: time(s) %lf\n",timeCounter[2]);
	printf("Runs Test: time(s) %lf\n",timeCounter[3]);
	printf("LongestRunOfOnes Test: time(s) %lf\n",timeCounter[4]);
	printf("Rank Test: time(s) %lf\n",timeCounter[5]);
	printf("DiscreteFourierTransform Test: time(s) %lf\n",timeCounter[6]);	
	printf("NonOverlappingTemplateMatchings Test: time(s) %lf\n",timeCounter[7]);
	printf("OverlappingTemplateMatchings Test: time(s) %lf\n",timeCounter[8]);
	printf("Universal Test: time(s) %lf\n",timeCounter[9]);
	printf("ApproximateEntropy Test: time(s) %lf\n",timeCounter[10]);
	printf("RandomExcursions Test: time(s) %lf\n",timeCounter[11]);
	printf("RandomExcursionsVariant Test: time(s) %lf\n",timeCounter[12]);
	printf("Serial Test: time(s) %lf\n",timeCounter[13]);
	printf("LinearComplexity Test: time(s) %lf\n",timeCounter[14]);
	if(testDM){
		printf("DeviationMeasureTotal Test: time(s) %lf\n",timeCounter[15]);
	}else{
		printf("DeviationMeasure1 Test: time(s) %lf\n",timeCounter[15]);
		printf("DeviationMeasure2 Test: time(s) %lf\n",timeCounter[16]);
		printf("DeviationMeasure1CS Test: time(s) %lf\n",timeCounter[17]);
		printf("DeviationMeasure2CS Test: time(s) %lf\n",timeCounter[18]);				
	}
	if(testJumpHop){
		printf("JumpHopComplexity Test: time(s) %lf\n",timeCounter[19]);
	}else{
		printf("Oddhopsum Test: time(s) %lf\n",timeCounter[19]);
		printf("Evenhopsum Test: time(s) %lf\n",timeCounter[20]);
		printf("Jump Test: time(s) %lf\n",timeCounter[21]);
	}
	printf("SelfCorrelation Test: time(s) %lf\n",timeCounter[22]);
	printf("BinaryDerivate Test: time(s) %lf\n",timeCounter[23]);
	printf("PokerDetect Test: time(s) %lf\n",timeCounter[24]);
	printf("RunsDistribution Test: time(s) %lf\n",timeCounter[25]);

	printf("total: time(s) %lf\n\n",sum);
}
