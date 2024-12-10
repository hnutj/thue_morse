#ifndef _INTERFACE_H_
#define _INTERFACE_H_

extern "C"{
    void NonOverlappingHelper(int m, int n,int templateNum, unsigned char* byteEpsilon, unsigned short* templates, unsigned int*container);
    void cufft(int n, double* input, double* output);
}

#endif