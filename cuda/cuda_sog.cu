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

#define SIZE_3GRAM_TABLE 0x1000000
#define CHAR_WIDTH_3GRAM 8

#define GET3GRAM(address) ((((uint32_t) (address)[0])) + (((uint32_t)((address)[1])) << CHAR_WIDTH_3GRAM) + (((uint32_t)((address)[2])) << (CHAR_WIDTH_3GRAM << 1)))

#define GET32(address) (((uint32_t)((address)[0]) << 24) + ((uint32_t)((address)[1]) << 16) + ((uint32_t)((address)[2]) << 8) + (address)[3])

texture<uint8_t, cudaTextureType1D> tex_T8;
texture<uint32_t, cudaTextureType1D> tex_scanner_hs;
texture<int, cudaTextureType1D> tex_scanner_index;
texture<uint8_t, cudaTextureType1D> tex_scanner_hs2;

//texture<unsigned char, cudaTextureType2D> tex_pattern;

__device__ int tex_sog_rkbt_verification8 ( unsigned char *d_pattern, size_t pattern_pitch, unsigned char *text, int pos, int m, int p_size ) {

	const uint8_t mask[] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};
	
	int i;

	uint32_t hs = GET32((text + pos)) ^ GET32((text + pos + 4));
	uint16_t hs2level = (uint16_t) ((hs >> 16) ^ hs);
	
	// check 2-level hash
	if ( tex1Dfetch ( tex_scanner_hs2, hs2level >> 3 ) & mask[hs2level & 0x07] ) {

		int lo = 0;
		int hi = p_size - 1;
		int mid;
		uint32_t hs_pat;
			
		// do the binary search
		while ( hi >= lo ) {

			mid = ( lo + hi ) / 2;
			hs_pat = tex1Dfetch( tex_scanner_hs, mid );

			if ( hs > hs_pat )
				lo = ++mid;

			else if ( hs < hs_pat )
				hi = --mid;
			
			//if text hash equals pattern hash verify the match
			else {
				// check for duplicates and patterns with same hash
				while ( mid > 0 && hs == tex1Dfetch( tex_scanner_hs, mid - 1 ) )
					mid--;
				
				do {
					for ( i = 0; i < m; i++ )
						//if ( tex2D ( tex_pattern, i, tex1Dfetch ( tex_scanner_index, mid ) ) != text[pos + i] )
						if ( d_pattern[tex1Dfetch( tex_scanner_index, mid ) * pattern_pitch + i] != text[pos + i] )
							break;
					
					if ( i == m )
						return 1;

					mid++;
				
				} while ( mid < p_size && hs == tex1Dfetch( tex_scanner_hs, mid ) );
			
				break;
			}
		}
	}
	return -1;
}

__device__ int sog_rkbt_verification8 ( uint32_t *scanner_hs, int *scanner_index, uint8_t *scanner_hs2, unsigned char *pattern, size_t pattern_pitch, unsigned char *text, int pos, int m, int p_size ) {

	const uint8_t mask[] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};
	
	int i;

	uint32_t hs = GET32((text + pos)) ^ GET32((text + pos + 4));
	uint16_t hs2level = (uint16_t) ((hs >> 16) ^ hs);
	
	// check 2-level hash
	if ( scanner_hs2[hs2level >> 3] & mask[hs2level & 0x07] ) {

		int lo = 0;
		int hi = p_size - 1;
		int mid;
		uint32_t hs_pat;
			
		// do the binary search
		while ( hi >= lo ) {

			mid = ( lo + hi ) / 2;
			hs_pat = scanner_hs[mid];

			if ( hs > hs_pat )
				lo = ++mid;

			else if ( hs < hs_pat )
				hi = --mid;
			
			//if text hash equals pattern hash verify the match
			else {
				// check for duplicates and patterns with same hash
				while ( mid > 0 && hs == scanner_hs[mid - 1] )
					mid--;
				
				do {
					for ( i = 0; i < m; i++ )
						if ( pattern[scanner_index[mid] * pattern_pitch + i] != text[pos + i] )
							break;
					
					if ( i == m )
						return 1;

					mid++;
				
				} while ( mid < p_size && hs == scanner_hs[mid] );
			
				break;
			}
		}
	}
	return -1;
}

__global__ void sog_kernel5 ( unsigned char *d_pattern, unsigned char *d_text, unsigned int *d_out, size_t pattern_pitch, int m, int n, int p_size, int numBlocks, int B, int sharedMemSize ) {
	
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
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
	
	int i, j, column, matches = 0;
	
	register uint8_t E = 0xff;
	
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
		
		for ( column = startThread; column < stopThread && column < sharedMemSize - 2; column++ ) {
		
			E = (E << 1) | tex1Dfetch ( tex_T8, GET3GRAM( s_array + column ) );
		
			if ( E & 0x20 )
				continue;
		
			if ( tex_sog_rkbt_verification8 ( d_pattern, pattern_pitch, s_array, column - m + B, m, p_size ) != -1 )
				matches++;
		}
		
		__syncthreads();
	}
	
	d_out[idx] = matches;
}

extern "C" void cuda_sog5 ( uint8_t *T8, uint32_t *scanner_hs, int *scanner_index, uint8_t *scanner_hs2, unsigned char *pattern, int m, unsigned char *text, int n, int p_size, int B ) {

	//Pointer for device memory
	unsigned char *d_pattern, *d_text;
	
	unsigned int *d_out;
	
	int *d_scanner_index;
	
	uint8_t *d_T8, *d_scanner_hs2;
	
	uint32_t *d_scanner_hs;
	
	size_t pattern_pitch;
	
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
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_hs, p_size * sizeof ( uint32_t ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_index, p_size * sizeof ( int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_hs2, 256 * 32 * sizeof ( uint8_t ) ) );
	
	//Allocate 2D device memory
	checkCudaErrors ( cudaMallocPitch ( &d_pattern, &pattern_pitch, m * sizeof ( unsigned char ), p_size ) );
	
	//Copy 1D host memory to device
	checkCudaErrors ( cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_T8, T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_hs, scanner_hs, p_size * sizeof ( uint32_t ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_index, scanner_index, p_size * sizeof ( int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_hs2, scanner_hs2, 256 * 32 * sizeof ( uint8_t ), cudaMemcpyHostToDevice ) );

	//Copy 2D host memory to device
	checkCudaErrors ( cudaMemcpy2D ( d_pattern, pattern_pitch, pattern, m * sizeof ( unsigned char ), m * sizeof ( unsigned char ), p_size, cudaMemcpyHostToDevice ) );
	
	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc int_desc = cudaCreateChannelDesc<int>();
	//cudaChannelFormatDesc char_desc = cudaCreateChannelDesc<unsigned char>();
	cudaChannelFormatDesc uint8_desc = cudaCreateChannelDesc<uint8_t>();
	
	checkCudaErrors ( cudaBindTexture ( 0, tex_T8, d_T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_hs, d_scanner_hs, p_size * sizeof ( uint32_t ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_index, d_scanner_index, p_size * sizeof ( int ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_hs2, d_scanner_hs2, 256 * 32 * sizeof ( uint8_t ) ) );
	
	//checkCudaErrors ( cudaBindTexture2D ( 0, tex_pattern, d_pattern, char_desc, m, p_size, pattern_pitch ) );
	
	//Create timer
	cudaEvent_t start, stop;

	float time;
	
	//Create the timer events
	cudaEventCreate ( &start );
	cudaEventCreate ( &stop );
	
	//Start the event clock	
	cudaEventRecord ( start, 0 );
	
	//Executing kernel in the device
	sog_kernel5<<<dimGrid, dimBlock, sharedMemSize + 16 * ( ( m - 1 ) / 16 + 1 )>>>( d_pattern, d_text, d_out, pattern_pitch, m, n, p_size, numBlocks, B, sharedMemSize );
	
	checkCUDAError("kernel invocation");

	cudaEventRecord ( stop, 0 );

	cudaEventSynchronize ( stop );
	
	cudaEventElapsedTime ( &time, start, stop );
	
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
	cudaFree ( d_pattern );
	cudaFree ( d_out );
	cudaFree ( d_T8 );
	cudaFree ( d_scanner_hs );
	cudaFree ( d_scanner_index );
	cudaFree ( d_scanner_hs2 );
}

__global__ void sog_kernel4 ( unsigned char *d_pattern, unsigned char *d_text, unsigned int *d_out, size_t pattern_pitch, int m, int n, int p_size, int numBlocks, int B, int sharedMemSize ) {
	
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
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
	
	int i, j, column;
	
	register uint8_t E = 0xff;
	
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
		
		for ( column = startThread; column < stopThread && column < sharedMemSize - 2; column++ ) {
		
			E = (E << 1) | tex1Dfetch ( tex_T8, GET3GRAM( s_array + column ) );
		
			if ( E & 0x20 )
				continue;
		
			if ( tex_sog_rkbt_verification8 ( d_pattern, pattern_pitch, s_array, column - m + B, m, p_size ) != -1 )
				d_out[idx]++;
		}
		
		__syncthreads();
	}
}

extern "C" void cuda_sog4 ( uint8_t *T8, uint32_t *scanner_hs, int *scanner_index, uint8_t *scanner_hs2, unsigned char *pattern, int m, unsigned char *text, int n, int p_size, int B ) {

	//Pointer for device memory
	unsigned char *d_pattern, *d_text;
	
	unsigned int *d_out;
	
	int *d_scanner_index;
	
	uint8_t *d_T8, *d_scanner_hs2;
	
	uint32_t *d_scanner_hs;
	
	size_t pattern_pitch;
	
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
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_hs, p_size * sizeof ( uint32_t ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_index, p_size * sizeof ( int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_hs2, 256 * 32 * sizeof ( uint8_t ) ) );
	
	//Allocate 2D device memory
	checkCudaErrors ( cudaMallocPitch ( &d_pattern, &pattern_pitch, m * sizeof ( unsigned char ), p_size ) );
	
	//Copy 1D host memory to device
	checkCudaErrors ( cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_T8, T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_hs, scanner_hs, p_size * sizeof ( uint32_t ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_index, scanner_index, p_size * sizeof ( int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_hs2, scanner_hs2, 256 * 32 * sizeof ( uint8_t ), cudaMemcpyHostToDevice ) );

	//Copy 2D host memory to device
	checkCudaErrors ( cudaMemcpy2D ( d_pattern, pattern_pitch, pattern, m * sizeof ( unsigned char ), m * sizeof ( unsigned char ), p_size, cudaMemcpyHostToDevice ) );
	
	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc int_desc = cudaCreateChannelDesc<int>();
	//cudaChannelFormatDesc char_desc = cudaCreateChannelDesc<unsigned char>();
	cudaChannelFormatDesc uint8_desc = cudaCreateChannelDesc<uint8_t>();
	
	checkCudaErrors ( cudaBindTexture ( 0, tex_T8, d_T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_hs, d_scanner_hs, p_size * sizeof ( uint32_t ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_index, d_scanner_index, p_size * sizeof ( int ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_hs2, d_scanner_hs2, 256 * 32 * sizeof ( uint8_t ) ) );
	
	//checkCudaErrors ( cudaBindTexture2D ( 0, tex_pattern, d_pattern, char_desc, m, p_size, pattern_pitch ) );
	
	//Create timer
	cudaEvent_t start, stop;

	float time;
	
	//Create the timer events
	cudaEventCreate ( &start );
	cudaEventCreate ( &stop );
	
	//Start the event clock	
	cudaEventRecord ( start, 0 );
	
	//Executing kernel in the device
	sog_kernel4<<<dimGrid, dimBlock, sharedMemSize + 16 * ( ( m - 1 ) / 16 + 1 )>>>( d_pattern, d_text, d_out, pattern_pitch, m, n, p_size, numBlocks, B, sharedMemSize );
	
	checkCUDAError("kernel invocation");

	cudaEventRecord ( stop, 0 );

	cudaEventSynchronize ( stop );
	
	cudaEventElapsedTime ( &time, start, stop );
	
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
	cudaFree ( d_pattern );
	cudaFree ( d_out );
	cudaFree ( d_T8 );
	cudaFree ( d_scanner_hs );
	cudaFree ( d_scanner_index );
	cudaFree ( d_scanner_hs2 );
}

__global__ void sog_kernel3 ( unsigned char *d_pattern, unsigned char *d_text, unsigned int *d_out, size_t pattern_pitch, int m, int n, int p_size, int numBlocks, int B, int sharedMemSize ) {
	
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	int charactersPerThread = sharedMemSize / blockDim.x;
	
	int startThread = charactersPerThread * threadIdx.x;
	int stopThread = startThread + charactersPerThread + m - B - 2;

	//Define space in shared memory
	extern __shared__ unsigned char s_array[];
	
	int i, j, column;
	
	register uint8_t E = 0xff;
	
	for ( int globalMemIndex = blockIdx.x * sharedMemSize; globalMemIndex < n; globalMemIndex += numBlocks * sharedMemSize ) {
		
		for ( i = globalMemIndex + threadIdx.x, j = 0 + threadIdx.x; ( j < sharedMemSize + m - 1 && i < n ); i+=blockDim.x, j+=blockDim.x )
			s_array[j] = d_text[i];
			
		__syncthreads();
		
		for ( column = startThread; column < stopThread && column < sharedMemSize - 2; column++ ) {
		
			E = (E << 1) | tex1Dfetch ( tex_T8, GET3GRAM( s_array + column ) );
		
			if ( E & 0x20 )
				continue;
		
			if ( tex_sog_rkbt_verification8 ( d_pattern, pattern_pitch, s_array, column - m + B, m, p_size ) != -1 )
				d_out[idx]++;
		}
		
		__syncthreads();
	}
}

extern "C" void cuda_sog3 ( uint8_t *T8, uint32_t *scanner_hs, int *scanner_index, uint8_t *scanner_hs2, unsigned char *pattern, int m, unsigned char *text, int n, int p_size, int B ) {

	//Pointer for device memory
	unsigned char *d_pattern, *d_text;
	
	unsigned int *d_out;
	
	int *d_scanner_index;
	
	uint8_t *d_T8, *d_scanner_hs2;
	
	uint32_t *d_scanner_hs;
	
	size_t pattern_pitch;
	
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
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_hs, p_size * sizeof ( uint32_t ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_index, p_size * sizeof ( int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_hs2, 256 * 32 * sizeof ( uint8_t ) ) );
	
	//Allocate 2D device memory
	checkCudaErrors ( cudaMallocPitch ( &d_pattern, &pattern_pitch, m * sizeof ( unsigned char ), p_size ) );
	
	//Copy 1D host memory to device
	checkCudaErrors ( cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_T8, T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_hs, scanner_hs, p_size * sizeof ( uint32_t ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_index, scanner_index, p_size * sizeof ( int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_hs2, scanner_hs2, 256 * 32 * sizeof ( uint8_t ), cudaMemcpyHostToDevice ) );

	//Copy 2D host memory to device
	checkCudaErrors ( cudaMemcpy2D ( d_pattern, pattern_pitch, pattern, m * sizeof ( unsigned char ), m * sizeof ( unsigned char ), p_size, cudaMemcpyHostToDevice ) );
	
	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc int_desc = cudaCreateChannelDesc<int>();
	//cudaChannelFormatDesc char_desc = cudaCreateChannelDesc<unsigned char>();
	cudaChannelFormatDesc uint8_desc = cudaCreateChannelDesc<uint8_t>();
	
	checkCudaErrors ( cudaBindTexture ( 0, tex_T8, d_T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_hs, d_scanner_hs, p_size * sizeof ( uint32_t ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_index, d_scanner_index, p_size * sizeof ( int ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_hs2, d_scanner_hs2, 256 * 32 * sizeof ( uint8_t ) ) );
	
	//checkCudaErrors ( cudaBindTexture2D ( 0, tex_pattern, d_pattern, char_desc, m, p_size, pattern_pitch ) );
	
	//Create timer
	cudaEvent_t start, stop;

	float time;
	
	//Create the timer events
	cudaEventCreate ( &start );
	cudaEventCreate ( &stop );
	
	//Start the event clock	
	cudaEventRecord ( start, 0 );
	
	//Executing kernel in the device
	sog_kernel3<<<dimGrid, dimBlock, sharedMemSize + m - 1>>>( d_pattern, d_text, d_out, pattern_pitch, m, n, p_size, numBlocks, B, sharedMemSize );
	
	checkCUDAError("kernel invocation");

	cudaEventRecord ( stop, 0 );

	cudaEventSynchronize ( stop );
	
	cudaEventElapsedTime ( &time, start, stop );
	
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
	cudaFree ( d_pattern );
	cudaFree ( d_out );
	cudaFree ( d_T8 );
	cudaFree ( d_scanner_hs );
	cudaFree ( d_scanner_index );
	cudaFree ( d_scanner_hs2 );
}

__global__ void sog_kernel2 ( unsigned char *d_pattern, unsigned char *d_text, unsigned int *d_out, size_t pattern_pitch, int m, int n, int p_size, int numBlocks, int B ) {
	
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	int charactersPerBlock = n / numBlocks;
	
	int startBlock = blockIdx.x * charactersPerBlock;
	int stopBlock = startBlock + charactersPerBlock;
	
	int charactersPerThread = ( stopBlock - startBlock ) / blockDim.x;
	
	int startThread = startBlock + charactersPerThread * threadIdx.x;
	int stopThread = startThread + charactersPerThread - B + m;
	
	int column;
	
	register uint8_t E = 0xff;
	
	for ( column = startThread; column < stopThread; column++ ) {
	
		if ( column < n - 2 ) {
		
			E = (E << 1) | tex1Dfetch ( tex_T8, GET3GRAM( d_text + column ) );
		
			if ( E & 0x20 )
				continue;
		
			if ( tex_sog_rkbt_verification8 ( d_pattern, pattern_pitch, d_text, column - m + B, m, p_size ) != -1 )
				d_out[idx]++;
		}
	}
}

extern "C" void cuda_sog2 ( uint8_t *T8, uint32_t *scanner_hs, int *scanner_index, uint8_t *scanner_hs2, unsigned char *pattern, int m, unsigned char *text, int n, int p_size, int B ) {

	//Pointer for device memory
	unsigned char *d_pattern, *d_text;
	
	unsigned int *d_out;
	
	int *d_scanner_index;
	
	uint8_t *d_T8, *d_scanner_hs2;
	
	uint32_t *d_scanner_hs;
	
	size_t pattern_pitch;
	
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
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_hs, p_size * sizeof ( uint32_t ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_index, p_size * sizeof ( int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_hs2, 256 * 32 * sizeof ( uint8_t ) ) );
	
	//Allocate 2D device memory
	checkCudaErrors ( cudaMallocPitch ( &d_pattern, &pattern_pitch, m * sizeof ( unsigned char ), p_size ) );
	
	//Copy 1D host memory to device
	checkCudaErrors ( cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_T8, T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_hs, scanner_hs, p_size * sizeof ( uint32_t ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_index, scanner_index, p_size * sizeof ( int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_hs2, scanner_hs2, 256 * 32 * sizeof ( uint8_t ), cudaMemcpyHostToDevice ) );

	//Copy 2D host memory to device
	checkCudaErrors ( cudaMemcpy2D ( d_pattern, pattern_pitch, pattern, m * sizeof ( unsigned char ), m * sizeof ( unsigned char ), p_size, cudaMemcpyHostToDevice ) );
	
	//Bind the preprocessing tables to the texture cache
	//cudaChannelFormatDesc int_desc = cudaCreateChannelDesc<int>();
	//cudaChannelFormatDesc char_desc = cudaCreateChannelDesc<unsigned char>();
	//cudaChannelFormatDesc uint8_desc = cudaCreateChannelDesc<uint8_t>();
	
	checkCudaErrors ( cudaBindTexture ( 0, tex_T8, d_T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_hs, d_scanner_hs, p_size * sizeof ( uint32_t ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_index, d_scanner_index, p_size * sizeof ( int ) ) );
	checkCudaErrors ( cudaBindTexture ( 0, tex_scanner_hs2, d_scanner_hs2, 256 * 32 * sizeof ( uint8_t ) ) );
	
	//checkCudaErrors ( cudaBindTexture2D ( 0, tex_pattern, d_pattern, char_desc, m, p_size, pattern_pitch ) );
	
	//Create timer
	cudaEvent_t start, stop;

	float time;
	
	//Create the timer events
	cudaEventCreate ( &start );
	cudaEventCreate ( &stop );
	
	//Start the event clock	
	cudaEventRecord ( start, 0 );
	
	//Executing kernel in the device
	sog_kernel2<<<dimGrid, dimBlock>>>( d_pattern, d_text, d_out, pattern_pitch, m, n, p_size, numBlocks, B );
	
	checkCUDAError("kernel invocation");

	cudaEventRecord ( stop, 0 );

	cudaEventSynchronize ( stop );
	
	cudaEventElapsedTime ( &time, start, stop );
	
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
	cudaFree ( d_pattern );
	cudaFree ( d_out );
	cudaFree ( d_T8 );
	cudaFree ( d_scanner_hs );
	cudaFree ( d_scanner_index );
	cudaFree ( d_scanner_hs2 );
}

__global__ void sog_kernel1 ( uint8_t *d_T8, uint32_t *d_scanner_hs, int *d_scanner_index, uint8_t *d_scanner_hs2, unsigned char *d_pattern, unsigned char *d_text, unsigned int *d_out, size_t pattern_pitch, int m, int n, int p_size, int numBlocks, int B ) {
	
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	int charactersPerBlock = n / numBlocks;
	
	int startBlock = blockIdx.x * charactersPerBlock;
	int stopBlock = startBlock + charactersPerBlock;
	
	int charactersPerThread = ( stopBlock - startBlock ) / blockDim.x;
	
	int startThread = startBlock + charactersPerThread * threadIdx.x;
	int stopThread = startThread + charactersPerThread - B + m; //-B + 1 + m - 1
	
	int column;
	
	register uint8_t E = 0xff;
	
	for ( column = startThread; column < stopThread; column++ ) {
	
		if ( column < n - 2 ) {
		
			E = (E << 1) | d_T8[GET3GRAM( d_text + column )];
		
			if ( E & 0x20 )
				continue;
		
			if ( sog_rkbt_verification8 ( d_scanner_hs, d_scanner_index, d_scanner_hs2, d_pattern, pattern_pitch, d_text, column - m + B, m, p_size ) != -1 )
				d_out[idx]++;
		}
	}
}
	
extern "C" void cuda_sog1 ( uint8_t *T8, uint32_t *scanner_hs, int *scanner_index, uint8_t *scanner_hs2, unsigned char *pattern, int m, unsigned char *text, int n, int p_size, int B ) {

	//Pointer for device memory
	unsigned char *d_pattern, *d_text;
	
	unsigned int *d_out;
	
	int *d_scanner_index;
	
	uint8_t *d_T8, *d_scanner_hs2;
	
	uint32_t *d_scanner_hs;
	
	size_t pattern_pitch;
	
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
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_hs, p_size * sizeof ( uint32_t ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_index, p_size * sizeof ( int ) ) );
	checkCudaErrors ( cudaMalloc ( ( void** ) &d_scanner_hs2, 256 * 32 * sizeof ( uint8_t ) ) );
	
	//Allocate 2D device memory
	checkCudaErrors ( cudaMallocPitch ( &d_pattern, &pattern_pitch, m * sizeof ( unsigned char ), p_size ) );
	
	//Copy 1D host memory to device
	checkCudaErrors ( cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_T8, T8, SIZE_3GRAM_TABLE * sizeof ( uint8_t ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_hs, scanner_hs, p_size * sizeof ( uint32_t ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_index, scanner_index, p_size * sizeof ( int ), cudaMemcpyHostToDevice ) );
	checkCudaErrors ( cudaMemcpy ( d_scanner_hs2, scanner_hs2, 256 * 32 * sizeof ( uint8_t ), cudaMemcpyHostToDevice ) );

	//Copy 2D host memory to device
	checkCudaErrors ( cudaMemcpy2D ( d_pattern, pattern_pitch, pattern, m * sizeof ( unsigned char ), m * sizeof ( unsigned char ), p_size, cudaMemcpyHostToDevice ) );
	
	//Create timer
	cudaEvent_t start, stop;

	float time;
	
	//Create the timer events
	cudaEventCreate ( &start );
	cudaEventCreate ( &stop );
	
	//Start the event clock	
	cudaEventRecord ( start, 0 );
	
	//Executing kernel in the device
	sog_kernel1<<<dimGrid, dimBlock>>>( d_T8, d_scanner_hs, d_scanner_index, d_scanner_hs2, d_pattern, d_text, d_out, pattern_pitch, m, n, p_size, numBlocks, B );
	
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
  	
	printf ("Kernel 1 matches \t%i\t time \t%f\n", matches, time/1000);
		
	//Free host and device memory
	free ( h_out );
	
	cudaFree ( d_text );
	cudaFree ( d_pattern );
	cudaFree ( d_out );
	cudaFree ( d_T8 );
	cudaFree ( d_scanner_hs );
	cudaFree ( d_scanner_index );
	cudaFree ( d_scanner_hs2 );
}

