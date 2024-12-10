#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <device_launch_parameters.h>
#include "../include/interface.h"
#include "../include/common.cuh"

cudaError_t ErrorCheck(cudaError_t error_code, const char* filename, int lineNumber)
{
    if (error_code != cudaSuccess)
    {
        printf("CUDA error:\r\ncode=%d, name=%s, description=%s\r\nfile=%s, line%d\r\n",
                error_code, cudaGetErrorName(error_code), cudaGetErrorString(error_code), filename, lineNumber);
        return error_code;
    }
    return error_code;
}

void cufft(int n, double* input, double* output){

    cudaSetDevice(0);
    int i;
    cufftDoubleReal *idata_cpu;
    cufftDoubleComplex *odata_cpu;
    cufftDoubleReal *idata_device;
    cufftDoubleComplex *odata_device;
    
    idata_cpu=(cufftDoubleReal*)malloc(n*sizeof(cufftDoubleReal));
    odata_cpu=(cufftDoubleComplex*)malloc((n/2+1)*sizeof(cufftDoubleComplex));
    for(i=0;i<n;i++)
    {
        idata_cpu[i]=input[i];
    }
    
    // host to device
    // ErrorCheck(cudaFree(0),__FILE__, __LINE__);
    ErrorCheck(cudaMalloc((void**)&idata_device,n*sizeof(cufftDoubleReal)),__FILE__, __LINE__);
    ErrorCheck(cudaMalloc((void**)&odata_device,(n/2+1)*sizeof(cufftDoubleComplex)),__FILE__, __LINE__);

    // 使用流化还可以再提高效率
    ErrorCheck(cudaMemcpy(idata_device,idata_cpu,n*sizeof(cufftDoubleReal),cudaMemcpyHostToDevice),__FILE__, __LINE__);

    // exec fft     
    // 单精度使用接口cufftExecR2C，双精度为D2Z
    cufftHandle plan;
    cufftPlan1d(&plan,n,CUFFT_D2Z,1);
    cufftExecD2Z(plan,(cufftDoubleReal*)idata_device,(cufftDoubleComplex*)odata_device);
    cudaDeviceSynchronize();
    
    // device to host
    ErrorCheck(cudaMemcpy(odata_cpu,odata_device,(n/2+1)*sizeof(cufftDoubleComplex),cudaMemcpyDeviceToHost),__FILE__, __LINE__);

    // for(i=0;i<n/2+1;i++)
    // {
    //     printf("%lf",odata_cpu[i].x);
    //     if(odata_cpu[i].y != 0.0 )
    //     {
    //         printf("+%lfi",odata_cpu[i].y);
    //     } 
    //     printf("\n");
    // }

	output[0] = sqrt(odata_cpu[0].x * odata_cpu[0].x);	    /* COMPUTE MAGNITUDE */
	
    //这边也可以用并行，速度会更快
	for ( i=1; i<n/2+1; i++ )
		output[i] = sqrt(pow(odata_cpu[i].x,2)+pow(odata_cpu[i].y,2)); 

    // for(i=0;i<n/2+1;i++)
    // {
    //     printf("%lf ",output[i]);
    // }
    // printf("\n");

    cufftDestroy(plan);
    
    free(idata_cpu);
    free(odata_cpu);
    cudaFree(odata_device);
    cudaFree(idata_device);
}