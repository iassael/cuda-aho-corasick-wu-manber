#ifndef CUDA_H
#define CUDA_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <cuda.h>

#define MAX(a,b) (a>b)?a:b

static void checkCUDAError(const char *msg) {

	cudaError_t err = cudaGetLastError();

	if (cudaSuccess != err) {
		fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
}

// This will output the proper CUDA error strings in the event that a CUDA host call returns an error
#define checkCudaErrors(err)           __checkCudaErrors (err, __FILE__, __LINE__)

inline static void __checkCudaErrors(cudaError err, const char *file,
		const int line) {

	if (cudaSuccess != err) {
		fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n", file, line,
				(int) err, cudaGetErrorString(err));
		exit(-1);
	}
}

#endif
