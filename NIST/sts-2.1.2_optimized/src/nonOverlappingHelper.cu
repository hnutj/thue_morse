#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/interface.h"
#include "../include/common.cuh"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          N O N O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// 手工实现(htonl(((int*)(ByteEpsilon+i))[0]
__device__ unsigned int myHtonls(unsigned int num1, unsigned int num2, int numOff1) 
{
	unsigned int result = 0;
	int numOff2 = numOff1 - 1;
	num1 >>= 8 * numOff1;
	for (int i = numOff1; i <= 3; i++) {
		result = result << 8 | (num1 & 0b11111111);
		num1 >>= 8;
	}
	for (int i = numOff2; i >= 0; i--) {
		result = result << 8 | (num2 & 0b11111111);
		num2 >>= 8;
	}
	return result;
}

//核函数
__global__ void NonOverlapping(int m, int n,unsigned char *cudaByteEpsilon,unsigned short *cudeTemplates, unsigned int *cudaContainer)
{
    const int jj = threadIdx.x;
    const int j = blockIdx.x;

	int N = 8;
	int M = n/N;
    unsigned short mask = powf(2, m) - 1;
    unsigned short B = cudeTemplates[jj]; //置于内层循环内外均可,置于内层为的是折叠循环便于omp并行
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
                test = (cudaByteEpsilon[p] >> k) & mask;
            }
            else if (k >= -8) {
                //涉及2bytes
                test = ((cudaByteEpsilon[p] << (-k)) | (cudaByteEpsilon[p + 1] >> (8 + k))) & mask;
            }
            else if (k >= -15) {
                //涉及3bytes
				int numOff1 = p % 4;
				int num1 = ((int*)(cudaByteEpsilon + p / 4 * 4))[0];
				int num2 = ((int*)(cudaByteEpsilon + (p / 4 + 1) * 4))[0];
				test = (myHtonls(num1, num2, numOff1) >> (24 + k)) & mask;			
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
    cudaContainer[jj*N+j] = W_obs;    
}

void NonOverlappingHelper(int m, int n, int templateNum, unsigned char* byteEpsilon, unsigned short* templates, unsigned int*container)
{
	// 锁定分块数为8
	int N = 8;
	// 每个子块的长度
	// int M = n/N;

    // 1、设置GPU设备
    // setGPU();
    cudaSetDevice(0);
    
    // 2、分配主机内存和设备内存，并初始化
    int byteNum=n / 8;
    int byteTotalNum = n % 8 == 0 ? byteNum + 4 : byteNum + 5;
    int stBytesCount1 = byteTotalNum * sizeof(unsigned char);
    int stBytesCount2 = templateNum * sizeof(unsigned short);
    int stBytesCount3 = templateNum * N * sizeof(unsigned int);

    unsigned char *cudaByteEpsilon;
    unsigned short *cudeTemplates;
    unsigned int *cudaContainer;

    //一次初始化后，程序cudamalloc()分配的内存不释放，继续使用，所有程序运行结束后，再一起释放
    ErrorCheck(cudaFree(0),__FILE__, __LINE__);

    ErrorCheck(cudaMalloc((void**)&cudaByteEpsilon, stBytesCount1),__FILE__, __LINE__);
    ErrorCheck(cudaMalloc((void**)&cudeTemplates, stBytesCount2),__FILE__, __LINE__);
    ErrorCheck(cudaMalloc((void**)&cudaContainer, stBytesCount3),__FILE__, __LINE__);

    // if (cudaByteEpsilon != NULL && cudeTemplates != NULL && cudaContainer != NULL)
    // {
    //     ErrorCheck(cudaMemset(cudaByteEpsilon, 0,stBytesCount1),__FILE__, __LINE__);  // 设备内存初始化为0
    //     ErrorCheck(cudaMemset(cudeTemplates, 0, stBytesCount2),__FILE__, __LINE__);
    //     ErrorCheck(cudaMemset(cudaContainer, 0, stBytesCount3),__FILE__, __LINE__);
    // }
    // else
    // {
    //     printf("fail to allocate memory\n");
    //     free(cudaByteEpsilon);
    //     free(cudeTemplates);
    //     free(cudaContainer);
    //     exit(-1);
    // }

    // 3、数据从主机复制到设备
    ErrorCheck(cudaMemcpy(cudaByteEpsilon, byteEpsilon, stBytesCount1, cudaMemcpyHostToDevice),__FILE__, __LINE__); 
    ErrorCheck(cudaMemcpy(cudeTemplates, templates, stBytesCount2, cudaMemcpyHostToDevice),__FILE__, __LINE__); 

    // 4、调用核函数
    dim3 block(templateNum);
    dim3 grid(N);
    NonOverlapping<<<grid, block>>>(m, n, cudaByteEpsilon, cudeTemplates, cudaContainer);   
    ErrorCheck(cudaGetLastError(), __FILE__, __LINE__);
    ErrorCheck(cudaDeviceSynchronize(), __FILE__, __LINE__);

    // 5、将计算得到的数据从设备传给主机
    ErrorCheck(cudaMemcpy(container, cudaContainer, stBytesCount3, cudaMemcpyDeviceToHost),__FILE__, __LINE__);

    // 6、释放内存
    // ErrorCheck(cudaFree(cudaByteEpsilon),__FILE__, __LINE__);
    // ErrorCheck(cudaFree(cudeTemplates),__FILE__, __LINE__);
    // ErrorCheck(cudaFree(cudaContainer),__FILE__, __LINE__);

    // ErrorCheck(cudaDeviceReset(),__FILE__, __LINE__);	
}
