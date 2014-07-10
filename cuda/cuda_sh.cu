/*This file is part of "A Hybrid Parallel Implementation of the Aho-Corasick and Wu-Manber Algorithms Using NVIDIA CUDA and MPI Evaluated on a Biological Sequence Database".

"A Hybrid Parallel Implementation of the Aho-Corasick and Wu-Manber Algorithms Using NVIDIA CUDA and MPI Evaluated on a Biological Sequence Database" is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"A Hybrid Parallel Implementation of the Aho-Corasick and Wu-Manber Algorithms Using NVIDIA CUDA and MPI Evaluated on a Biological Sequence Database" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with "A Hybrid Parallel Implementation of the Aho-Corasick and Wu-Manber Algorithms Using NVIDIA CUDA and MPI Evaluated on a Biological Sequence Database".  If not, see <http://www.gnu.org/licenses/>.*/

#include "cuda.h"

texture<int, cudaTextureType2D> tex_state_transition;
texture<unsigned int, cudaTextureType1D> tex_state_final;
texture<int, cudaTextureType1D> tex_bmBc;

//Optimization 4: Group the write to global memory operations
__global__ void sh_kernel5 ( unsigned char *d_text, unsigned int *d_out, int m, int n, int p_size, int alphabet, int numBlocks, int sharedMemSize ) {

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int r, s;
	
	int i, j, column, matches = 0;
	
	int charactersPerThread = sharedMemSize / blockDim.x;
	
	int startThread = charactersPerThread * threadIdx.x + m - 1;
	int stopThread = startThread + charactersPerThread;
	
	//Define space in shared memory
	extern __shared__ unsigned char s_array[];
	
	//cast data to uint4
	uint4 *uint4_text = reinterpret_cast < uint4 * > ( d_text );
	uint4 uint4_var;
	
	//recast data to uchar4
	uchar4 c0, c4, c8, c12;
	
	for ( int globalMemIndex = blockIdx.x * sharedMemSize; globalMemIndex < n; globalMemIndex += numBlocks * sharedMemSize ) {
	
		for ( i = globalMemIndex/16 + threadIdx.x, j = 0 + threadIdx.x; j < sharedMemSize / 16 && i < n / 16; i+=blockDim.x, j+=blockDim.x ) {
			
			uint4_var = uint4_text[i];
			
			//recast data back to char after the memory transaction
			c0 = *reinterpret_cast<uchar4 *> ( &uint4_var.x );
			c4 = *reinterpret_cast<uchar4 *> ( &uint4_var.y );
			c8 = *reinterpret_cast<uchar4 *> ( &uint4_var.z );
			c12 = *reinterpret_cast<uchar4 *> ( &uint4_var.w );

			s_array[j * 16 + 0] = c0.x;
                        s_array[j * 16 + 1] = c0.y;
                        s_array[j * 16 + 2] = c0.z;
                        s_array[j * 16 + 3] = c0.w;
                        
                        s_array[j * 16 + 4] = c4.x;
                        s_array[j * 16 + 5] = c4.y;
                        s_array[j * 16 + 6] = c4.z;
                        s_array[j * 16 + 7] = c4.w;
                        
                        s_array[j * 16 + 8] = c8.x;
                        s_array[j * 16 + 9] = c8.y;
                        s_array[j * 16 + 10] = c8.z;
                        s_array[j * 16 + 11] = c8.w;
                        
                        s_array[j * 16 + 12] = c12.x;
                        s_array[j * 16 + 13] = c12.y;
                        s_array[j * 16 + 14] = c12.z;
                        s_array[j * 16 + 15] = c12.w;
		}

		//Add m - 1 redundant characters at the end of the shared memory
		//FIXME: optimize this!!
		if ( threadIdx.x < m - 1 )
			s_array[sharedMemSize + threadIdx.x] = d_text[globalMemIndex + sharedMemSize + threadIdx.x];
			
		__syncthreads();
		
		column = startThread;
	
		while ( column < stopThread && globalMemIndex + column < n  ) {
	
			r = 0;
			j = 0;

			while ( j < m && ( s = tex2D ( tex_state_transition, s_array[column - j], r ) ) != -1 ) {
				
				r = s;
				j++;
			}
		
			matches += tex1Dfetch ( tex_state_final, r );
		
			column += tex1Dfetch ( tex_bmBc, d_text[column] );
		}
		
		__syncthreads();
	}
	
	d_out[idx] = matches;
}

extern "C" void cuda_sh5 ( int m, unsigned char *text, int n, int p_size, int alphabet, int *state_transition, unsigned int *state_final, int *bmBc ) {

	//Pointer for device memory
	int *d_state_transition, *d_bmBc;
	unsigned int *d_state_final, *d_out;
	
	unsigned char *d_text;

	size_t pitch;
	
	int numBlocks = 30, numThreadsPerBlock = 256, sharedMemSize = 16128;
	dim3 dimGrid ( numBlocks );
	dim3 dimBlock ( numThreadsPerBlock );
	
	if ( n < numBlocks * numThreadsPerBlock * m ) {
		printf("The text size is too small\n");
		exit(1);
	}
	
	//Allocate host memory for results array
	unsigned int *h_out = ( unsigned int * ) malloc ( numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) );
	memset ( h_out, 0, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) );
	
	//Allocate 1D device memory
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_text, n * sizeof ( unsigned char ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_bmBc, alphabet * sizeof ( int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ) );
	
	//Allocate 2D device memory
	checkCudaErrors ( cudaMallocPitch ( &d_state_transition, &pitch, alphabet * sizeof ( int ), ( m * p_size + 1 ) ) );
	
	//Copy 1D host memory to device
	checkCudaErrors ( cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_state_final, state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_bmBc, bmBc, alphabet * sizeof ( int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	
	//Copy 2D host memory to device
	checkCudaErrors ( cudaMemcpy2D ( d_state_transition, pitch, state_transition, alphabet * sizeof ( int ), alphabet * sizeof ( int ), ( m * p_size + 1 ), cudaMemcpyHostToDevice ) );
	
	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc desc = cudaCreateChannelDesc<int>();
	checkCudaErrors ( cudaBindTexture2D ( 0, tex_state_transition, d_state_transition, desc, alphabet, m * p_size + 1, pitch ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_state_final, d_state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_bmBc, d_bmBc, alphabet * sizeof ( int ) ) );
	
	//Create timer
	cudaEvent_t start, stop;

	float time;

	//Create the timer events
	cudaEventCreate ( &start );
	cudaEventCreate ( &stop );
	
	//Start the event clock	
	cudaEventRecord ( start, 0 );
	
	//Executing kernel in the device
	sh_kernel5<<<dimGrid, dimBlock, sharedMemSize + 16 * ( ( m - 1 ) / 16 + 1 )>>>( d_text, d_out, m, n, p_size, alphabet, numBlocks, sharedMemSize );
	checkCUDAError("kernel invocation");
	
	cudaEventRecord ( stop, 0 );

	cudaEventSynchronize ( stop );
	
	cudaEventElapsedTime ( &time, start, stop );
	
	cudaEventDestroy ( start );
	cudaEventDestroy ( stop );

	//Get back the results from the device
	cudaMemcpy ( h_out, d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyDeviceToHost );
	   
  	//Look at the results
  	int i, matches = 0;
  	
  	for ( i = 0; i < numBlocks * numThreadsPerBlock; i++ )
  		matches += h_out[i];
  	
	printf ("Kernel 5 matches \t%i\t time \t%f\n", matches, time/1000);
	
	//Free host and device memory
	free ( h_out );

	cudaFree ( d_text );
	cudaFree ( d_state_transition );
	cudaFree ( d_state_final );
	cudaFree ( d_bmBc );	
	cudaFree ( d_out );
}

//Optimization 3: Read 16 byte words per thread with coalescing. Uint4 words are extracted to shared memory after fetching from global memory so it is not as efficient as extracting on a per thread basis as with the ac uint4 optimization.
__global__ void sh_kernel4 ( unsigned char *d_text, unsigned int *d_out, int m, int n, int p_size, int alphabet, int numBlocks, int sharedMemSize ) {

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int r, s;
	
	int i, j, column;
	
	int charactersPerThread = sharedMemSize / blockDim.x;
	
	int startThread = charactersPerThread * threadIdx.x + m - 1;
	int stopThread = startThread + charactersPerThread;
	
	//Define space in shared memory
	extern __shared__ unsigned char s_array[];
	
	//cast data to uint4
	uint4 *uint4_text = reinterpret_cast < uint4 * > ( d_text );
	uint4 uint4_var;
	
	//recast data to uchar4
	uchar4 c0, c4, c8, c12;
	
	for ( int globalMemIndex = blockIdx.x * sharedMemSize; globalMemIndex < n; globalMemIndex += numBlocks * sharedMemSize ) {
	
		for ( i = globalMemIndex/16 + threadIdx.x, j = 0 + threadIdx.x; j < sharedMemSize / 16 && i < n / 16; i+=blockDim.x, j+=blockDim.x ) {
			
			uint4_var = uint4_text[i];
			
			//recast data back to char after the memory transaction
			c0 = *reinterpret_cast<uchar4 *> ( &uint4_var.x );
			c4 = *reinterpret_cast<uchar4 *> ( &uint4_var.y );
			c8 = *reinterpret_cast<uchar4 *> ( &uint4_var.z );
			c12 = *reinterpret_cast<uchar4 *> ( &uint4_var.w );

			s_array[j * 16 + 0] = c0.x;
                        s_array[j * 16 + 1] = c0.y;
                        s_array[j * 16 + 2] = c0.z;
                        s_array[j * 16 + 3] = c0.w;
                        
                        s_array[j * 16 + 4] = c4.x;
                        s_array[j * 16 + 5] = c4.y;
                        s_array[j * 16 + 6] = c4.z;
                        s_array[j * 16 + 7] = c4.w;
                        
                        s_array[j * 16 + 8] = c8.x;
                        s_array[j * 16 + 9] = c8.y;
                        s_array[j * 16 + 10] = c8.z;
                        s_array[j * 16 + 11] = c8.w;
                        
                        s_array[j * 16 + 12] = c12.x;
                        s_array[j * 16 + 13] = c12.y;
                        s_array[j * 16 + 14] = c12.z;
                        s_array[j * 16 + 15] = c12.w;
		}

		//Add m - 1 redundant characters at the end of the shared memory
		//FIXME: optimize this!!
		if ( threadIdx.x < m - 1 )
			s_array[sharedMemSize + threadIdx.x] = d_text[globalMemIndex + sharedMemSize + threadIdx.x];
			
		__syncthreads();
		
		column = startThread;
	
		while ( column < stopThread && globalMemIndex + column < n  ) {
	
			r = 0;
			j = 0;

			while ( j < m && ( s = tex2D ( tex_state_transition, s_array[column - j], r ) ) != -1 ) {
				
				r = s;
				j++;
			}
		
			d_out[idx] += tex1Dfetch ( tex_state_final, r );
		
			column += tex1Dfetch ( tex_bmBc, d_text[column] );
		}
		
		__syncthreads();
	}
}

extern "C" void cuda_sh4 ( int m, unsigned char *text, int n, int p_size, int alphabet, int *state_transition, unsigned int *state_final, int *bmBc ) {

	//Pointer for device memory
	int *d_state_transition, *d_bmBc;
	unsigned int *d_state_final, *d_out;
	
	unsigned char *d_text;

	size_t pitch;
	
	int numBlocks = 30, numThreadsPerBlock = 256, sharedMemSize = 16128;
	dim3 dimGrid ( numBlocks );
	dim3 dimBlock ( numThreadsPerBlock );
	
	if ( n < numBlocks * numThreadsPerBlock * m ) {
		printf("The text size is too small\n");
		exit(1);
	}
	
	//Allocate host memory for results array
	unsigned int *h_out = ( unsigned int * ) malloc ( numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) );
	memset ( h_out, 0, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) );
	
	//Allocate 1D device memory
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_text, n * sizeof ( unsigned char ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_bmBc, alphabet * sizeof ( int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ) );
	
	//Allocate 2D device memory
	checkCudaErrors ( cudaMallocPitch ( &d_state_transition, &pitch, alphabet * sizeof ( int ), ( m * p_size + 1 ) ) );
	
	//Copy 1D host memory to device
	checkCudaErrors ( cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_state_final, state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_bmBc, bmBc, alphabet * sizeof ( int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	
	//Copy 2D host memory to device
	checkCudaErrors ( cudaMemcpy2D ( d_state_transition, pitch, state_transition, alphabet * sizeof ( int ), alphabet * sizeof ( int ), ( m * p_size + 1 ), cudaMemcpyHostToDevice ) );
	
	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc desc = cudaCreateChannelDesc<int>();
	checkCudaErrors ( cudaBindTexture2D ( 0, tex_state_transition, d_state_transition, desc, alphabet, m * p_size + 1, pitch ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_state_final, d_state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_bmBc, d_bmBc, alphabet * sizeof ( int ) ) );
	
	//Create timer
	cudaEvent_t start, stop;

	float time;

	//Create the timer events
	cudaEventCreate ( &start );
	cudaEventCreate ( &stop );
	
	//Start the event clock	
	cudaEventRecord ( start, 0 );
	
	//Executing kernel in the device
	sh_kernel4<<<dimGrid, dimBlock, sharedMemSize + 16 * ( ( m - 1 ) / 16 + 1 )>>>( d_text, d_out, m, n, p_size, alphabet, numBlocks, sharedMemSize );
	checkCUDAError("kernel invocation");
	
	cudaEventRecord ( stop, 0 );

	cudaEventSynchronize ( stop );
	
	cudaEventElapsedTime ( &time, start, stop );
	
	cudaEventDestroy ( start );
	cudaEventDestroy ( stop );

	//Get back the results from the device
	cudaMemcpy ( h_out, d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyDeviceToHost );
	   
  	//Look at the results
  	int i, matches = 0;
  	
  	for ( i = 0; i < numBlocks * numThreadsPerBlock; i++ )
  		matches += h_out[i];
  	
	printf ("Kernel 4 matches \t%i\t time \t%f\n", matches, time/1000);
	
	//Free host and device memory
	free ( h_out );

	cudaFree ( d_text );
	cudaFree ( d_state_transition );
	cudaFree ( d_state_final );
	cudaFree ( d_bmBc );	
	cudaFree ( d_out );
}

//Optimization 2: Read sharedMemSize characters byte-to-byte from global memory to shared memory to coalescelce memory transactions 
__global__ void sh_kernel3 ( unsigned char *d_text, unsigned int *d_out, int m, int n, int p_size, int alphabet, int numBlocks, int sharedMemSize ) {

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int r, s;
	
	int i, j, column;
	
	int charactersPerThread = sharedMemSize / blockDim.x;
	
	int startThread = charactersPerThread * threadIdx.x + m - 1;
	int stopThread = startThread + charactersPerThread;
	
	//Define space in shared memory
	extern __shared__ unsigned char s_array[];
	
	for ( int globalMemIndex = blockIdx.x * sharedMemSize; globalMemIndex < n; globalMemIndex += numBlocks * sharedMemSize ) {
	
	
		for ( i = globalMemIndex + threadIdx.x, j = 0 + threadIdx.x; ( j < sharedMemSize + m - 1 && i < n ); i+=blockDim.x, j+=blockDim.x )
			s_array[j] = d_text[i];
			
		__syncthreads();
		
		column = startThread;
	
		while ( column < stopThread && globalMemIndex + column < n  ) {
	
			r = 0;
			j = 0;

			while ( j < m && ( s = tex2D ( tex_state_transition, s_array[column - j], r ) ) != -1 ) {
				
				r = s;
				j++;
			}
		
			d_out[idx] += tex1Dfetch ( tex_state_final, r );
		
			column += tex1Dfetch ( tex_bmBc, d_text[column] );
		}
		
		__syncthreads();
	}
}

extern "C" void cuda_sh3 ( int m, unsigned char *text, int n, int p_size, int alphabet, int *state_transition, unsigned int *state_final, int *bmBc ) {

	//Pointer for device memory
	int *d_state_transition, *d_bmBc;
	unsigned int *d_state_final, *d_out;
	
	unsigned char *d_text;

	size_t pitch;
	
	int numBlocks = 30, numThreadsPerBlock = 256, sharedMemSize = 16128;
	dim3 dimGrid ( numBlocks );
	dim3 dimBlock ( numThreadsPerBlock );
	
	if ( n < numBlocks * numThreadsPerBlock * m ) {
		printf("The text size is too small\n");
		exit(1);
	}
	
	//Allocate host memory for results array
	unsigned int *h_out = ( unsigned int * ) malloc ( numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) );
	memset ( h_out, 0, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) );
	
	//Allocate 1D device memory
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_text, n * sizeof ( unsigned char ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_bmBc, alphabet * sizeof ( int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ) );
	
	//Allocate 2D device memory
	checkCudaErrors ( cudaMallocPitch ( &d_state_transition, &pitch, alphabet * sizeof ( int ), ( m * p_size + 1 ) ) );
	
	//Copy 1D host memory to device
	checkCudaErrors ( cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_state_final, state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_bmBc, bmBc, alphabet * sizeof ( int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	
	//Copy 2D host memory to device
	checkCudaErrors ( cudaMemcpy2D ( d_state_transition, pitch, state_transition, alphabet * sizeof ( int ), alphabet * sizeof ( int ), ( m * p_size + 1 ), cudaMemcpyHostToDevice ) );
	
	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc desc = cudaCreateChannelDesc<int>();
	checkCudaErrors ( cudaBindTexture2D ( 0, tex_state_transition, d_state_transition, desc, alphabet, m * p_size + 1, pitch ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_state_final, d_state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_bmBc, d_bmBc, alphabet * sizeof ( int ) ) );
	
	//Create timer
	cudaEvent_t start, stop;

	float time;

	//Create the timer events
	cudaEventCreate ( &start );
	cudaEventCreate ( &stop );
	
	//Start the event clock	
	cudaEventRecord ( start, 0 );
	
	//Executing kernel in the device
	sh_kernel3<<<dimGrid, dimBlock, sharedMemSize + m - 1>>>( d_text, d_out, m, n, p_size, alphabet, numBlocks, sharedMemSize );
	checkCUDAError("kernel invocation");
	
	cudaEventRecord ( stop, 0 );

	cudaEventSynchronize ( stop );
	
	cudaEventElapsedTime ( &time, start, stop );
	
	cudaEventDestroy ( start );
	cudaEventDestroy ( stop );

	//Get back the results from the device
	cudaMemcpy ( h_out, d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyDeviceToHost );
	   
  	//Look at the results
  	int i, matches = 0;
  	
  	for ( i = 0; i < numBlocks * numThreadsPerBlock; i++ )
  		matches += h_out[i];
  	
	printf ("Kernel 3 matches \t%i\t time \t%f\n", matches, time/1000);
	
	//Free host and device memory
	free ( h_out );

	cudaFree ( d_text );
	cudaFree ( d_state_transition );
	cudaFree ( d_state_final );
	cudaFree ( d_bmBc );
	cudaFree ( d_out );
}

//Optimization 1: Use the texture cache for the pattern
__global__ void sh_kernel2 ( unsigned char *d_text, unsigned int *d_out, int m, int n, int p_size, int alphabet, int numBlocks ) {

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	int charactersPerBlock = n / numBlocks;
	
	int startBlock = blockIdx.x * charactersPerBlock;
	int stopBlock = startBlock + charactersPerBlock;
	
	int charactersPerThread = ( stopBlock - startBlock ) / blockDim.x;
	
	int startThread = startBlock + charactersPerThread * threadIdx.x + m - 1;
	int stopThread = startThread + charactersPerThread;

	int r, s;
	
	int column = startThread, j;
	
	while ( column < stopThread ) {
	
		r = 0;
		j = 0;

		while ( j < m && ( s = tex2D ( tex_state_transition, d_text[column - j], r ) ) != -1 ) {
				
			r = s;
			j++;
		}
		
		d_out[idx] += tex1Dfetch ( tex_state_final, r );
		
		column += tex1Dfetch ( tex_bmBc, d_text[column] );
	}
}

extern "C" void cuda_sh2 ( int m, unsigned char *text, int n, int p_size, int alphabet, int *state_transition, unsigned int *state_final, int *bmBc ) {

	//Pointer for device memory
	int *d_state_transition, *d_bmBc;
	unsigned int *d_state_final, *d_out;
	
	unsigned char *d_text;
	
	size_t pitch;
	
	int numBlocks = 30, numThreadsPerBlock = 256;
	dim3 dimGrid ( numBlocks );
	dim3 dimBlock ( numThreadsPerBlock );
	
	if ( n < numBlocks * numThreadsPerBlock * m ) {
		printf("The text size is too small\n");
		exit(1);
	}
	
	//Allocate host memory for results array
	unsigned int *h_out = ( unsigned int * ) malloc ( numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) );
	memset ( h_out, 0, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) );

	//Allocate 1D device memory
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_text, n * sizeof ( unsigned char ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_bmBc, alphabet * sizeof ( int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ) );
	
	//Allocate 2D device memory
	checkCudaErrors ( cudaMallocPitch ( &d_state_transition, &pitch, alphabet * sizeof ( int ), ( m * p_size + 1 ) ) );
	
	//Copy 1D host memory to device
	checkCudaErrors ( cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_state_final, state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_bmBc, bmBc, alphabet * sizeof ( int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );

	//Copy 2D host memory to device
	checkCudaErrors ( cudaMemcpy2D ( d_state_transition, pitch, state_transition, alphabet * sizeof ( int ), alphabet * sizeof ( int ), ( m * p_size + 1 ), cudaMemcpyHostToDevice ) );
	
	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc desc = cudaCreateChannelDesc<int>();
	checkCudaErrors ( cudaBindTexture2D ( 0, tex_state_transition, d_state_transition, desc, alphabet, m * p_size + 1, pitch ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_state_final, d_state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_bmBc, d_bmBc, alphabet * sizeof ( int ) ) );
	
	//Create timer
	cudaEvent_t start, stop;

	float time;

	//Create the timer events
	cudaEventCreate ( &start );
	cudaEventCreate ( &stop );
	
	//Start the event clock	
	cudaEventRecord ( start, 0 );
	
	//Executing kernel in the device
	sh_kernel2<<<dimGrid, dimBlock>>>( d_text, d_out, m, n, p_size, alphabet, numBlocks );
	
	checkCUDAError("kernel invocation");
	
	cudaEventRecord ( stop, 0 );

	cudaEventSynchronize ( stop );
	
	cudaEventElapsedTime ( &time, start, stop );
	
	cudaEventDestroy ( start );
	cudaEventDestroy ( stop );

	//Get back the results from the device
	cudaMemcpy ( h_out, d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyDeviceToHost );
   
  	//Look at the results
  	int i, matches = 0;
  	for ( i = 0; i < numBlocks * numThreadsPerBlock; i++ )
  		matches += h_out[i];
  	
	printf ("Kernel 2 matches \t%i\t time \t%f\n", matches, time/1000);
		
	//Free host and device memory
	free ( h_out );
	
	cudaFree ( d_text );
	cudaFree ( d_state_transition );
	cudaFree ( d_state_final );
	cudaFree ( d_bmBc );
	cudaFree ( d_out );
}

__global__ void sh_kernel1 ( int *d_state_transition, unsigned int *d_state_final, int *d_bmBc, unsigned char *d_text, unsigned int *d_out, size_t pitch, int m, int n, int p_size, int alphabet, int numBlocks ) {
	
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int effective_pitch = pitch / sizeof ( int );
	
	int charactersPerBlock = n / numBlocks;
	
	int startBlock = blockIdx.x * charactersPerBlock;
	int stopBlock = startBlock + charactersPerBlock;
	
	int charactersPerThread = ( stopBlock - startBlock ) / blockDim.x;
	
	int startThread = startBlock + charactersPerThread * threadIdx.x + m - 1;
	int stopThread = startThread + charactersPerThread;

	int r, s;
	
	int column = startThread, j;
	
	while ( column < stopThread ) {
	
		r = 0;
		j = 0;

		while ( j < m && ( s = d_state_transition[r * effective_pitch + d_text[column - j]] ) != -1 ) {
				
			r = s;
			j++;
		}
		
		d_out[idx] += d_state_final[r];
		
		column += d_bmBc[d_text[column]];
	}
}
	
extern "C" void cuda_sh1 ( int m, unsigned char *text, int n, int p_size, int alphabet, int *state_transition, unsigned int *state_final, int *bmBc ) {

	//Pointer for device memory
	int *d_state_transition, *d_bmBc;
	unsigned int *d_state_final, *d_out;
	
	unsigned char *d_text;

	size_t pitch;
	
	int numBlocks = 30, numThreadsPerBlock = 256;
	dim3 dimGrid ( numBlocks );
	dim3 dimBlock ( numThreadsPerBlock );
	
	if ( n < numBlocks * numThreadsPerBlock * m ) {
		printf("The text size is too small\n");
		exit(1);
	}
	
	//Allocate host memory for results array
	unsigned int *h_out = ( unsigned int * ) malloc ( numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) );
	memset ( h_out, 0, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) );

	//Allocate 1D device memory
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_text, n * sizeof ( unsigned char ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_bmBc, alphabet * sizeof ( int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ) );
	
	//Allocate 2D device memory
	checkCudaErrors ( cudaMallocPitch ( &d_state_transition, &pitch, alphabet * sizeof ( int ), ( m * p_size + 1 ) ) );
	
	//Copy 1D host memory to device
	checkCudaErrors ( cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_state_final, state_final, ( m * p_size + 1 ) * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_bmBc, bmBc, alphabet * sizeof ( int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );

	//Copy 2D host memory to device
	checkCudaErrors ( cudaMemcpy2D ( d_state_transition, pitch, state_transition, alphabet * sizeof ( int ), alphabet * sizeof ( int ), ( m * p_size + 1 ), cudaMemcpyHostToDevice ) );
	
	//Create timer
	cudaEvent_t start, stop;

	float time;
	
	//Create the timer events
	cudaEventCreate ( &start );
	cudaEventCreate ( &stop );
	
	//Start the event clock	
	cudaEventRecord ( start, 0 );
	
	//cudaPrintfInit();
	
	//Executing kernel in the device
	sh_kernel1<<<dimGrid, dimBlock>>>( d_state_transition, d_state_final, d_bmBc, d_text, d_out, pitch, m, n, p_size, alphabet, numBlocks );
	
	checkCUDAError("kernel invocation");
	
	//cudaPrintfDisplay(stdout, true);
	//cudaPrintfEnd();

	cudaEventRecord ( stop, 0 );

	cudaEventSynchronize ( stop );
	
	cudaEventElapsedTime ( &time, start, stop );
	
	cudaEventDestroy ( start );
	cudaEventDestroy ( stop );
	
	//Get back the results from the device
	cudaMemcpy ( h_out, d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyDeviceToHost );
	
  	//Look at the results
  	int i, matches = 0;
  	for ( i = 0; i < numBlocks * numThreadsPerBlock; i++ )
  		matches += h_out[i];
  	
	printf ("Kernel 1 matches \t%i\t time \t%f\n", matches, time/1000);
		
	//Free host and device memory
	free ( h_out );
	
	cudaFree ( d_text );
	cudaFree ( d_state_transition );
	cudaFree ( d_state_final );
	cudaFree ( d_bmBc );
	cudaFree ( d_out );
}
