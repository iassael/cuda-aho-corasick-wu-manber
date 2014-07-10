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

texture<int, cudaTextureType1D> tex_SHIFT;
texture<int, cudaTextureType1D> tex_PREFIX_size;

texture<int, cudaTextureType2D> tex_PREFIX_value;
texture<int, cudaTextureType2D> tex_PREFIX_index;
texture<unsigned char, cudaTextureType2D> tex_pattern;

unsigned int determine_shiftsize(int alphabet) {

	//the maximum size of the hash value of the B-size suffix of the patterns for the Wu-Manber algorithm
	if (alphabet == 2)
		return 22; // 1 << 2 + 1 << 2 + 1 + 1

	else if (alphabet == 4)
		return 64; // 3 << 2 + 3 << 2 + 3 + 1

	else if (alphabet == 8)
		return 148; // 7 << 2 + 7 << 2 + 7 + 1

	else if (alphabet == 20)
		return 400; // 19 << 2 + 19 << 2 + 19 + 1

	else if (alphabet == 128)
		return 2668; // 127 << 2 + 127 << 2 + 127 + 1

	else if (alphabet == 256)
		return 5356; //304 << 2 + 304 << 2 + 304 + 1

	else if (alphabet == 512)
		return 10732; //560 << 2 + 560 << 2 + 560 + 1

	else if (alphabet == 1024)
		return 21484; //1072 << 2 + 1072 << 2 + 1072 + 1

	else {
		printf("The alphabet size is not supported by wu-manber\n");
		exit(1);
	}
}

//FIXME: Test the implementation with fewer patterns so as the preprocessing tables can fit inside the texture cache

__global__ void wm_kernel5(int *d_PREFIX_index, unsigned char *d_text,
		unsigned int *d_out, size_t PREFIX_index_pitch, int m, int n,
		int p_size, int alphabet, int numBlocks, int sharedMemSize) {

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int effective_PREFIX_index_pitch = PREFIX_index_pitch / sizeof(int);

	int charactersPerThread = sharedMemSize / blockDim.x;

	int startThread = charactersPerThread * threadIdx.x + m - 1;
	int stopThread = startThread + charactersPerThread;

	//Define space in shared memory
	extern __shared__ unsigned char s_array[];

	//cast data to uint4
	uint4 *uint4_text = reinterpret_cast<uint4 *>(d_text);
	uint4 uint4_var;

	//recast data to uchar4
	uchar4 c0, c4, c8, c12;

	unsigned short m_nBitsInShift = 2;

	int column, i, j, l, matches = 0;

	unsigned int hash1, hash2;

	size_t shift;

	for (int globalMemIndex = blockIdx.x * sharedMemSize; globalMemIndex < n;
			globalMemIndex += numBlocks * sharedMemSize) {

		for (i = globalMemIndex / 16 + threadIdx.x, j = 0 + threadIdx.x;
				j < sharedMemSize / 16 && i < n / 16; i += blockDim.x, j +=
						blockDim.x) {

			uint4_var = uint4_text[i];

			//recast data back to char after the memory transaction
			c0 = *reinterpret_cast<uchar4 *>(&uint4_var.x);
			c4 = *reinterpret_cast<uchar4 *>(&uint4_var.y);
			c8 = *reinterpret_cast<uchar4 *>(&uint4_var.z);
			c12 = *reinterpret_cast<uchar4 *>(&uint4_var.w);

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
		if (threadIdx.x < m - 1)
			s_array[sharedMemSize + threadIdx.x] = d_text[globalMemIndex
					+ sharedMemSize + threadIdx.x];

		__syncthreads();

		column = startThread;

		while (column < stopThread) {

			hash1 = s_array[column - 2];
			hash1 <<= m_nBitsInShift;
			hash1 += s_array[column - 1];
			hash1 <<= m_nBitsInShift;
			hash1 += s_array[column];

			shift = tex1Dfetch(tex_SHIFT, hash1);

			if (shift == 0) {

				hash2 = s_array[column - m + 1];
				hash2 <<= m_nBitsInShift;
				hash2 += s_array[column - m + 2];

				//For every pattern with the same suffix as the text
				for (i = 0; i < tex1Dfetch(tex_PREFIX_size, hash1); i++) {

					//If the prefix of the pattern matches that of the text
					if (hash2 == tex2D(tex_PREFIX_value, i, hash1)) {

						//memcmp implementation
						for (l = 0; l < m; l++)
							if (tex2D(tex_pattern, l,
									d_PREFIX_index[hash1
											* effective_PREFIX_index_pitch + i])
									!= s_array[column - m + 1 + l])
								break;

						if (l == m) {
							matches++;
							break;
						}
					}
				}

				column++;
			} else
				column += shift;
		}
		__syncthreads();
	}

	d_out[idx] = matches;
}

extern "C" int cuda_wm5(unsigned char *pattern, int m, unsigned char *text,
		int n, int p_size, int alphabet, int B, int *SHIFT, int *PREFIX_value,
		int *PREFIX_index, int *PREFIX_size, double *gpuTime) {

	//Pointer for device memory
	int *d_SHIFT, *d_PREFIX_value, *d_PREFIX_index, *d_PREFIX_size;

	unsigned char *d_pattern, *d_text;

	unsigned int *d_out;

	size_t pattern_pitch, PREFIX_value_pitch, PREFIX_index_pitch;

	int numBlocks = 30, numThreadsPerBlock = 256, sharedMemSize = 16128;
	dim3 dimGrid(numBlocks);
	dim3 dimBlock(numThreadsPerBlock);

	if (n < numBlocks * numThreadsPerBlock * m) {
		printf("The text size is too small\n");
		exit(1);
	}

	unsigned int shiftsize = determine_shiftsize(alphabet);

	//Allocate host memory for results array
	unsigned int *h_out = (unsigned int *) malloc(
			numBlocks * numThreadsPerBlock * sizeof(unsigned int));
	memset(h_out, 0, numBlocks * numThreadsPerBlock * sizeof(unsigned int));

	//Allocate 1D device memory
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_text, n * sizeof ( unsigned char ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_SHIFT, shiftsize * sizeof ( int ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_PREFIX_size, shiftsize * sizeof ( int ) ));

	//Allocate 2D device memory
	checkCudaErrors(
			cudaMallocPitch ( &d_pattern, &pattern_pitch, m * sizeof ( unsigned char ), p_size ));
	checkCudaErrors(
			cudaMallocPitch ( &d_PREFIX_value, &PREFIX_value_pitch, p_size * sizeof ( int ), shiftsize ));
	checkCudaErrors(
			cudaMallocPitch ( &d_PREFIX_index, &PREFIX_index_pitch, p_size * sizeof ( int ), shiftsize ));

	//Copy 1D host memory to device
	checkCudaErrors(
			cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_SHIFT, SHIFT, shiftsize * sizeof ( int ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_PREFIX_size, PREFIX_size, shiftsize * sizeof ( int ), cudaMemcpyHostToDevice ));

	//Copy 2D host memory to device
	checkCudaErrors(
			cudaMemcpy2D ( d_pattern, pattern_pitch, pattern, m * sizeof ( unsigned char ), m * sizeof ( unsigned char ), p_size, cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy2D ( d_PREFIX_value, PREFIX_value_pitch, PREFIX_value, p_size * sizeof ( int ), p_size * sizeof ( int ), shiftsize, cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy2D ( d_PREFIX_index, PREFIX_index_pitch, PREFIX_index, p_size * sizeof ( int ), p_size * sizeof ( int ), shiftsize, cudaMemcpyHostToDevice ));

	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc int_desc = cudaCreateChannelDesc<int>();
	cudaChannelFormatDesc char_desc = cudaCreateChannelDesc<unsigned char>();
	checkCudaErrors(
			cudaBindTexture ( 0, tex_SHIFT, d_SHIFT, shiftsize * sizeof ( int ) ));
	checkCudaErrors(
			cudaBindTexture ( 0, tex_PREFIX_size, d_PREFIX_size, shiftsize * sizeof ( int ) ));

	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_PREFIX_value, d_PREFIX_value, int_desc, p_size, shiftsize, PREFIX_value_pitch ));
	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_PREFIX_index, d_PREFIX_index, int_desc, p_size, shiftsize, PREFIX_index_pitch ));
	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_pattern, d_pattern, char_desc, m, p_size, pattern_pitch ));

	//Create timer
	cudaEvent_t start, stop;

	float time;

	//Create the timer events
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//Start the event clock	
	cudaEventRecord(start, 0);

	//Executing kernel in the device
	wm_kernel5<<<dimGrid, dimBlock, sharedMemSize + 16 * ((m - 1) / 16 + 1)>>>(
			d_PREFIX_index, d_text, d_out, PREFIX_index_pitch, m, n, p_size,
			alphabet, numBlocks, sharedMemSize);

	checkCUDAError("kernel invocation");

	cudaEventRecord(stop, 0);

	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	//Get back the results from the device
	cudaMemcpy(h_out, d_out,
			numBlocks * numThreadsPerBlock * sizeof(unsigned int),
			cudaMemcpyDeviceToHost);

	//Look at the results
	int i, matches = 0;
	for (i = 0; i < numBlocks * numThreadsPerBlock; i++)
		matches += h_out[i];

	//printf("Kernel 5 matches \t%i\t time \t%f\n", matches, time / 1000);
	*gpuTime = time / 1000;

	//Free host and device memory
	free(h_out);

	cudaFree(d_text);
	cudaFree(d_pattern);
	cudaFree(d_SHIFT);
	cudaFree(d_PREFIX_value);
	cudaFree(d_PREFIX_index);
	cudaFree(d_PREFIX_size);
	cudaFree(d_out);
	return matches;
}

__global__ void wm_kernel4(int *d_PREFIX_index, unsigned char *d_text,
		unsigned int *d_out, size_t PREFIX_index_pitch, int m, int n,
		int p_size, int alphabet, int numBlocks, int sharedMemSize) {

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int effective_PREFIX_index_pitch = PREFIX_index_pitch / sizeof(int);

	int charactersPerThread = sharedMemSize / blockDim.x;

	int startThread = charactersPerThread * threadIdx.x + m - 1;
	int stopThread = startThread + charactersPerThread;

	//Define space in shared memory
	extern __shared__ unsigned char s_array[];

	//cast data to uint4
	uint4 *uint4_text = reinterpret_cast<uint4 *>(d_text);
	uint4 uint4_var;

	//recast data to uchar4
	uchar4 c0, c4, c8, c12;

	unsigned short m_nBitsInShift = 2;

	int column, i, j, l;

	unsigned int hash1, hash2;

	size_t shift;

	for (int globalMemIndex = blockIdx.x * sharedMemSize; globalMemIndex < n;
			globalMemIndex += numBlocks * sharedMemSize) {

		for (i = globalMemIndex / 16 + threadIdx.x, j = 0 + threadIdx.x;
				j < sharedMemSize / 16 && i < n / 16; i += blockDim.x, j +=
						blockDim.x) {

			uint4_var = uint4_text[i];

			//recast data back to char after the memory transaction
			c0 = *reinterpret_cast<uchar4 *>(&uint4_var.x);
			c4 = *reinterpret_cast<uchar4 *>(&uint4_var.y);
			c8 = *reinterpret_cast<uchar4 *>(&uint4_var.z);
			c12 = *reinterpret_cast<uchar4 *>(&uint4_var.w);

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
		if (threadIdx.x < m - 1)
			s_array[sharedMemSize + threadIdx.x] = d_text[globalMemIndex
					+ sharedMemSize + threadIdx.x];

		__syncthreads();

		column = startThread;

		while (column < stopThread) {

			hash1 = s_array[column - 2];
			hash1 <<= m_nBitsInShift;
			hash1 += s_array[column - 1];
			hash1 <<= m_nBitsInShift;
			hash1 += s_array[column];

			shift = tex1Dfetch(tex_SHIFT, hash1);

			if (shift == 0) {

				hash2 = s_array[column - m + 1];
				hash2 <<= m_nBitsInShift;
				hash2 += s_array[column - m + 2];

				//For every pattern with the same suffix as the text
				for (i = 0; i < tex1Dfetch(tex_PREFIX_size, hash1); i++) {

					//If the prefix of the pattern matches that of the text
					if (hash2 == tex2D(tex_PREFIX_value, i, hash1)) {

						//memcmp implementation
						for (l = 0; l < m; l++)
							if (tex2D(tex_pattern, l,
									d_PREFIX_index[hash1
											* effective_PREFIX_index_pitch + i])
									!= s_array[column - m + 1 + l])
								break;

						if (l == m) {
							d_out[idx]++;
							break;
						}
					}
				}

				column++;
			} else
				column += shift;
		}
		__syncthreads();
	}
}

extern "C" int cuda_wm4(unsigned char *pattern, int m, unsigned char *text,
		int n, int p_size, int alphabet, int B, int *SHIFT, int *PREFIX_value,
		int *PREFIX_index, int *PREFIX_size, double *gpuTime) {

	//Pointer for device memory
	int *d_SHIFT, *d_PREFIX_value, *d_PREFIX_index, *d_PREFIX_size;

	unsigned char *d_pattern, *d_text;

	unsigned int *d_out;

	size_t pattern_pitch, PREFIX_value_pitch, PREFIX_index_pitch;

	int numBlocks = 30, numThreadsPerBlock = 256, sharedMemSize = 16128;
	dim3 dimGrid(numBlocks);
	dim3 dimBlock(numThreadsPerBlock);

	if (n < numBlocks * numThreadsPerBlock * m) {
		printf("The text size is too small\n");
		exit(1);
	}

	unsigned int shiftsize = determine_shiftsize(alphabet);

	//Allocate host memory for results array
	unsigned int *h_out = (unsigned int *) malloc(
			numBlocks * numThreadsPerBlock * sizeof(unsigned int));
	memset(h_out, 0, numBlocks * numThreadsPerBlock * sizeof(unsigned int));

	//Allocate 1D device memory
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_text, n * sizeof ( unsigned char ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_SHIFT, shiftsize * sizeof ( int ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_PREFIX_size, shiftsize * sizeof ( int ) ));

	//Allocate 2D device memory
	checkCudaErrors(
			cudaMallocPitch ( &d_pattern, &pattern_pitch, m * sizeof ( unsigned char ), p_size ));
	checkCudaErrors(
			cudaMallocPitch ( &d_PREFIX_value, &PREFIX_value_pitch, p_size * sizeof ( int ), shiftsize ));
	checkCudaErrors(
			cudaMallocPitch ( &d_PREFIX_index, &PREFIX_index_pitch, p_size * sizeof ( int ), shiftsize ));

	//Copy 1D host memory to device
	checkCudaErrors(
			cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_SHIFT, SHIFT, shiftsize * sizeof ( int ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_PREFIX_size, PREFIX_size, shiftsize * sizeof ( int ), cudaMemcpyHostToDevice ));

	//Copy 2D host memory to device
	checkCudaErrors(
			cudaMemcpy2D ( d_pattern, pattern_pitch, pattern, m * sizeof ( unsigned char ), m * sizeof ( unsigned char ), p_size, cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy2D ( d_PREFIX_value, PREFIX_value_pitch, PREFIX_value, p_size * sizeof ( int ), p_size * sizeof ( int ), shiftsize, cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy2D ( d_PREFIX_index, PREFIX_index_pitch, PREFIX_index, p_size * sizeof ( int ), p_size * sizeof ( int ), shiftsize, cudaMemcpyHostToDevice ));

	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc int_desc = cudaCreateChannelDesc<int>();
	cudaChannelFormatDesc char_desc = cudaCreateChannelDesc<unsigned char>();
	checkCudaErrors(
			cudaBindTexture ( 0, tex_SHIFT, d_SHIFT, shiftsize * sizeof ( int ) ));
	checkCudaErrors(
			cudaBindTexture ( 0, tex_PREFIX_size, d_PREFIX_size, shiftsize * sizeof ( int ) ));

	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_PREFIX_value, d_PREFIX_value, int_desc, p_size, shiftsize, PREFIX_value_pitch ));
	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_PREFIX_index, d_PREFIX_index, int_desc, p_size, shiftsize, PREFIX_index_pitch ));
	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_pattern, d_pattern, char_desc, m, p_size, pattern_pitch ));

	//Create timer
	cudaEvent_t start, stop;

	float time;

	//Create the timer events
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//Start the event clock	
	cudaEventRecord(start, 0);

	//Executing kernel in the device
	wm_kernel4<<<dimGrid, dimBlock, sharedMemSize + 16 * ((m - 1) / 16 + 1)>>>(
			d_PREFIX_index, d_text, d_out, PREFIX_index_pitch, m, n, p_size,
			alphabet, numBlocks, sharedMemSize);

	checkCUDAError("kernel invocation");

	cudaEventRecord(stop, 0);

	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	//Get back the results from the device
	cudaMemcpy(h_out, d_out,
			numBlocks * numThreadsPerBlock * sizeof(unsigned int),
			cudaMemcpyDeviceToHost);

	//Look at the results
	int i, matches = 0;
	for (i = 0; i < numBlocks * numThreadsPerBlock; i++)
		matches += h_out[i];

	//printf("Kernel 4 matches \t%i\t time \t%f\n", matches, time / 1000);
	*gpuTime = time / 1000;

	//Free host and device memory
	free(h_out);

	cudaFree(d_text);
	cudaFree(d_pattern);
	cudaFree(d_SHIFT);
	cudaFree(d_PREFIX_value);
	cudaFree(d_PREFIX_index);
	cudaFree(d_PREFIX_size);
	cudaFree(d_out);
	return matches;
}

__global__ void wm_kernel3(int *d_PREFIX_index, unsigned char *d_text,
		unsigned int *d_out, size_t PREFIX_index_pitch, int m, int n,
		int p_size, int alphabet, int numBlocks, int sharedMemSize) {

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int effective_PREFIX_index_pitch = PREFIX_index_pitch / sizeof(int);

	int charactersPerThread = sharedMemSize / blockDim.x;

	int startThread = charactersPerThread * threadIdx.x + m - 1;
	int stopThread = startThread + charactersPerThread;

	//Define space in shared memory
	extern __shared__ unsigned char s_array[];

	unsigned short m_nBitsInShift = 2;

	int column, i, j, l;

	unsigned int hash1, hash2;

	size_t shift;

	for (int globalMemIndex = blockIdx.x * sharedMemSize; globalMemIndex < n;
			globalMemIndex += numBlocks * sharedMemSize) {

		for (i = globalMemIndex + threadIdx.x, j = 0 + threadIdx.x;
				(j < sharedMemSize + m - 1 && i < n); i += blockDim.x, j +=
						blockDim.x)
			s_array[j] = d_text[i];

		__syncthreads();

		column = startThread;

		while (column < stopThread) {

			hash1 = s_array[column - 2];
			hash1 <<= m_nBitsInShift;
			hash1 += s_array[column - 1];
			hash1 <<= m_nBitsInShift;
			hash1 += s_array[column];

			shift = tex1Dfetch(tex_SHIFT, hash1);

			if (shift == 0) {

				hash2 = s_array[column - m + 1];
				hash2 <<= m_nBitsInShift;
				hash2 += s_array[column - m + 2];

				//For every pattern with the same suffix as the text
				for (i = 0; i < tex1Dfetch(tex_PREFIX_size, hash1); i++) {

					//If the prefix of the pattern matches that of the text
					if (hash2 == tex2D(tex_PREFIX_value, i, hash1)) {

						//memcmp implementation
						for (l = 0; l < m; l++)
							if (tex2D(tex_pattern, l,
									d_PREFIX_index[hash1
											* effective_PREFIX_index_pitch + i])
									!= s_array[column - m + 1 + l])
								break;

						if (l == m) {
							d_out[idx]++;
							break;
						}
					}
				}

				column++;
			} else
				column += shift;
		}
		__syncthreads();
	}
}

extern "C" int cuda_wm3(unsigned char *pattern, int m, unsigned char *text,
		int n, int p_size, int alphabet, int B, int *SHIFT, int *PREFIX_value,
		int *PREFIX_index, int *PREFIX_size, double *gpuTime) {

	//Pointer for device memory
	int *d_SHIFT, *d_PREFIX_value, *d_PREFIX_index, *d_PREFIX_size;

	unsigned char *d_pattern, *d_text;

	unsigned int *d_out;

	size_t pattern_pitch, PREFIX_value_pitch, PREFIX_index_pitch;

	int numBlocks = 30, numThreadsPerBlock = 256, sharedMemSize = 16128;
	dim3 dimGrid(numBlocks);
	dim3 dimBlock(numThreadsPerBlock);

	if (n < numBlocks * numThreadsPerBlock * m) {
		printf("The text size is too small\n");
		exit(1);
	}

	unsigned int shiftsize = determine_shiftsize(alphabet);

	//Allocate host memory for results array
	unsigned int *h_out = (unsigned int *) malloc(
			numBlocks * numThreadsPerBlock * sizeof(unsigned int));
	memset(h_out, 0, numBlocks * numThreadsPerBlock * sizeof(unsigned int));

	//Allocate 1D device memory
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_text, n * sizeof ( unsigned char ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_SHIFT, shiftsize * sizeof ( int ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_PREFIX_size, shiftsize * sizeof ( int ) ));

	//Allocate 2D device memory
	checkCudaErrors(
			cudaMallocPitch ( &d_pattern, &pattern_pitch, m * sizeof ( unsigned char ), p_size ));
	checkCudaErrors(
			cudaMallocPitch ( &d_PREFIX_value, &PREFIX_value_pitch, p_size * sizeof ( int ), shiftsize ));
	checkCudaErrors(
			cudaMallocPitch ( &d_PREFIX_index, &PREFIX_index_pitch, p_size * sizeof ( int ), shiftsize ));

	//Copy 1D host memory to device
	checkCudaErrors(
			cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_SHIFT, SHIFT, shiftsize * sizeof ( int ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_PREFIX_size, PREFIX_size, shiftsize * sizeof ( int ), cudaMemcpyHostToDevice ));

	//Copy 2D host memory to device
	checkCudaErrors(
			cudaMemcpy2D ( d_pattern, pattern_pitch, pattern, m * sizeof ( unsigned char ), m * sizeof ( unsigned char ), p_size, cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy2D ( d_PREFIX_value, PREFIX_value_pitch, PREFIX_value, p_size * sizeof ( int ), p_size * sizeof ( int ), shiftsize, cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy2D ( d_PREFIX_index, PREFIX_index_pitch, PREFIX_index, p_size * sizeof ( int ), p_size * sizeof ( int ), shiftsize, cudaMemcpyHostToDevice ));

	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc int_desc = cudaCreateChannelDesc<int>();
	cudaChannelFormatDesc char_desc = cudaCreateChannelDesc<unsigned char>();
	checkCudaErrors(
			cudaBindTexture ( 0, tex_SHIFT, d_SHIFT, shiftsize * sizeof ( int ) ));
	checkCudaErrors(
			cudaBindTexture ( 0, tex_PREFIX_size, d_PREFIX_size, shiftsize * sizeof ( int ) ));

	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_PREFIX_value, d_PREFIX_value, int_desc, p_size, shiftsize, PREFIX_value_pitch ));
	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_PREFIX_index, d_PREFIX_index, int_desc, p_size, shiftsize, PREFIX_index_pitch ));
	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_pattern, d_pattern, char_desc, m, p_size, pattern_pitch ));

	//Create timer
	cudaEvent_t start, stop;

	float time;

	//Create the timer events
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//Start the event clock	
	cudaEventRecord(start, 0);

	//Executing kernel in the device
	wm_kernel3<<<dimGrid, dimBlock, sharedMemSize + m - 1>>>(d_PREFIX_index,
			d_text, d_out, PREFIX_index_pitch, m, n, p_size, alphabet,
			numBlocks, sharedMemSize);

	checkCUDAError("kernel invocation");

	cudaEventRecord(stop, 0);

	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	//Get back the results from the device
	cudaMemcpy(h_out, d_out,
			numBlocks * numThreadsPerBlock * sizeof(unsigned int),
			cudaMemcpyDeviceToHost);

	//Look at the results
	int i, matches = 0;
	for (i = 0; i < numBlocks * numThreadsPerBlock; i++)
		matches += h_out[i];

	//printf("Kernel 3 matches \t%i\t time \t%f\n", matches, time / 1000);
	*gpuTime = time / 1000;

	//Free host and device memory
	free(h_out);

	cudaFree(d_text);
	cudaFree(d_pattern);
	cudaFree(d_SHIFT);
	cudaFree(d_PREFIX_value);
	cudaFree(d_PREFIX_index);
	cudaFree(d_PREFIX_size);
	cudaFree(d_out);
	return matches;
}

__global__ void wm_kernel2(int *d_PREFIX_index, unsigned char *d_text,
		unsigned int *d_out, size_t PREFIX_index_pitch, int m, int n,
		int p_size, int alphabet, int numBlocks) {

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int effective_PREFIX_index_pitch = PREFIX_index_pitch / sizeof(int);

	int charactersPerBlock = n / numBlocks;

	int startBlock = blockIdx.x * charactersPerBlock;
	int stopBlock = startBlock + charactersPerBlock;

	int charactersPerThread = (stopBlock - startBlock) / blockDim.x;

	int startThread = startBlock + charactersPerThread * threadIdx.x + m - 1;
	int stopThread = startThread + charactersPerThread;

	unsigned short m_nBitsInShift = 2;

	int column = startThread, i, l;

	unsigned int hash1, hash2;

	size_t shift;

	while (column < stopThread) {

		hash1 = d_text[column - 2];
		hash1 <<= m_nBitsInShift;
		hash1 += d_text[column - 1];
		hash1 <<= m_nBitsInShift;
		hash1 += d_text[column];

		shift = tex1Dfetch(tex_SHIFT, hash1);

		if (shift == 0) {

			hash2 = d_text[column - m + 1];
			hash2 <<= m_nBitsInShift;
			hash2 += d_text[column - m + 2];

			//For every pattern with the same suffix as the text
			for (i = 0; i < tex1Dfetch(tex_PREFIX_size, hash1); i++) {

				//If the prefix of the pattern matches that of the text
				if (hash2 == tex2D(tex_PREFIX_value, i, hash1)) {

					//memcmp implementation
					for (l = 0; l < m; l++)
						if (tex2D(tex_pattern, l,
								d_PREFIX_index[hash1
										* effective_PREFIX_index_pitch + i])
								!= d_text[column - m + 1 + l])
							break;

					if (l == m) {
						d_out[idx]++;
						break;
					}
				}
			}

			column++;
		} else
			column += shift;
	}
}

extern "C" int cuda_wm2(unsigned char *pattern, int m, unsigned char *text,
		int n, int p_size, int alphabet, int B, int *SHIFT, int *PREFIX_value,
		int *PREFIX_index, int *PREFIX_size, double *gpuTime) {

	//Pointer for device memory
	int *d_SHIFT, *d_PREFIX_value, *d_PREFIX_index, *d_PREFIX_size;

	unsigned char *d_pattern, *d_text;

	unsigned int *d_out;

	size_t pattern_pitch, PREFIX_value_pitch, PREFIX_index_pitch;

	int numBlocks = 30, numThreadsPerBlock = 256;
	dim3 dimGrid(numBlocks);
	dim3 dimBlock(numThreadsPerBlock);

	if (n < numBlocks * numThreadsPerBlock * m) {
		printf("The text size is too small\n");
		exit(1);
	}

	unsigned int shiftsize = determine_shiftsize(alphabet);

	//Allocate host memory for results array
	unsigned int *h_out = (unsigned int *) malloc(
			numBlocks * numThreadsPerBlock * sizeof(unsigned int));
	memset(h_out, 0, numBlocks * numThreadsPerBlock * sizeof(unsigned int));

	//Allocate 1D device memory
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_text, n * sizeof ( unsigned char ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_SHIFT, shiftsize * sizeof ( int ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_PREFIX_size, shiftsize * sizeof ( int ) ));

	//Allocate 2D device memory
	checkCudaErrors(
			cudaMallocPitch ( &d_pattern, &pattern_pitch, m * sizeof ( unsigned char ), p_size ));
	checkCudaErrors(
			cudaMallocPitch ( &d_PREFIX_value, &PREFIX_value_pitch, p_size * sizeof ( int ), shiftsize ));
	checkCudaErrors(
			cudaMallocPitch ( &d_PREFIX_index, &PREFIX_index_pitch, p_size * sizeof ( int ), shiftsize ));

	//Copy 1D host memory to device
	checkCudaErrors(
			cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_SHIFT, SHIFT, shiftsize * sizeof ( int ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_PREFIX_size, PREFIX_size, shiftsize * sizeof ( int ), cudaMemcpyHostToDevice ));

	//Copy 2D host memory to device
	checkCudaErrors(
			cudaMemcpy2D ( d_pattern, pattern_pitch, pattern, m * sizeof ( unsigned char ), m * sizeof ( unsigned char ), p_size, cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy2D ( d_PREFIX_value, PREFIX_value_pitch, PREFIX_value, p_size * sizeof ( int ), p_size * sizeof ( int ), shiftsize, cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy2D ( d_PREFIX_index, PREFIX_index_pitch, PREFIX_index, p_size * sizeof ( int ), p_size * sizeof ( int ), shiftsize, cudaMemcpyHostToDevice ));

	//Bind the preprocessing tables to the texture cache
	cudaChannelFormatDesc int_desc = cudaCreateChannelDesc<int>();
	cudaChannelFormatDesc char_desc = cudaCreateChannelDesc<unsigned char>();
	checkCudaErrors(
			cudaBindTexture ( 0, tex_SHIFT, d_SHIFT, shiftsize * sizeof ( int ) ));
	checkCudaErrors(
			cudaBindTexture ( 0, tex_PREFIX_size, d_PREFIX_size, shiftsize * sizeof ( int ) ));

	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_PREFIX_value, d_PREFIX_value, int_desc, p_size, shiftsize, PREFIX_value_pitch ));
	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_PREFIX_index, d_PREFIX_index, int_desc, p_size, shiftsize, PREFIX_index_pitch ));
	checkCudaErrors(
			cudaBindTexture2D ( 0, tex_pattern, d_pattern, char_desc, m, p_size, pattern_pitch ));

	//Create timer
	cudaEvent_t start, stop;

	float time;

	//Create the timer events
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//Start the event clock	
	cudaEventRecord(start, 0);

	//Executing kernel in the device
	wm_kernel2<<<dimGrid, dimBlock>>>(d_PREFIX_index, d_text, d_out,
			PREFIX_index_pitch, m, n, p_size, alphabet, numBlocks);

	checkCUDAError("kernel invocation");

	cudaEventRecord(stop, 0);

	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	//Get back the results from the device
	cudaMemcpy(h_out, d_out,
			numBlocks * numThreadsPerBlock * sizeof(unsigned int),
			cudaMemcpyDeviceToHost);

	//Look at the results
	int i, matches = 0;
	for (i = 0; i < numBlocks * numThreadsPerBlock; i++)
		matches += h_out[i];

	//printf("Kernel 2 matches \t%i\t time \t%f\n", matches, time / 1000);
	*gpuTime = time / 1000;

	//Free host and device memory
	free(h_out);

	cudaFree(d_text);
	cudaFree(d_pattern);
	cudaFree(d_SHIFT);
	cudaFree(d_PREFIX_value);
	cudaFree(d_PREFIX_index);
	cudaFree(d_PREFIX_size);
	cudaFree(d_out);
	return matches;
}

__global__ void wm_kernel1(int *d_SHIFT, int *d_PREFIX_value,
		int *d_PREFIX_index, int *d_PREFIX_size, unsigned char *d_pattern,
		unsigned char *d_text, unsigned int *d_out, size_t pattern_pitch,
		size_t PREFIX_value_pitch, size_t PREFIX_index_pitch, int m, int n,
		int p_size, int alphabet, int numBlocks) {

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int effective_PREFIX_value_pitch = PREFIX_value_pitch / sizeof(int);
	int effective_PREFIX_index_pitch = PREFIX_index_pitch / sizeof(int);

	int charactersPerBlock = n / numBlocks;

	int startBlock = blockIdx.x * charactersPerBlock;
	int stopBlock = startBlock + charactersPerBlock;

	int charactersPerThread = (stopBlock - startBlock) / blockDim.x;

	int startThread = startBlock + charactersPerThread * threadIdx.x + m - 1;
	int stopThread = startThread + charactersPerThread;

	unsigned short m_nBitsInShift = 2;

	int column = startThread, i, l;

	unsigned int hash1, hash2;

	size_t shift;

	while (column < stopThread) {

		hash1 = d_text[column - 2];
		hash1 <<= m_nBitsInShift;
		hash1 += d_text[column - 1];
		hash1 <<= m_nBitsInShift;
		hash1 += d_text[column];

		shift = d_SHIFT[hash1];

		if (shift == 0) {

			hash2 = d_text[column - m + 1];
			hash2 <<= m_nBitsInShift;
			hash2 += d_text[column - m + 2];

			//For every pattern with the same suffix as the text
			for (i = 0; i < d_PREFIX_size[hash1]; i++) {

				//If the prefix of the pattern matches that of the text
				if (hash2
						== d_PREFIX_value[hash1 * effective_PREFIX_value_pitch
								+ i]) {

					//memcmp implementation
					for (l = 0; l < m; l++)
						if (d_pattern[d_PREFIX_index[hash1
								* effective_PREFIX_index_pitch + i]
								* pattern_pitch + l]
								!= d_text[column - m + 1 + l])
							break;

					if (l == m) {
						d_out[idx]++;
						break;
					}
				}
			}

			column++;
		} else
			column += shift;
	}
}

extern "C" int cuda_wm1(unsigned char *pattern, int m, unsigned char *text,
		int n, int p_size, int alphabet, int B, int *SHIFT, int *PREFIX_value,
		int *PREFIX_index, int *PREFIX_size, double *gpuTime) {

	//Pointer for device memory
	int *d_SHIFT, *d_PREFIX_value, *d_PREFIX_index, *d_PREFIX_size;

	unsigned char *d_pattern, *d_text;

	unsigned int *d_out;

	size_t pattern_pitch, PREFIX_value_pitch, PREFIX_index_pitch;

	int numBlocks = 30, numThreadsPerBlock = 256;
	dim3 dimGrid(numBlocks);
	dim3 dimBlock(numThreadsPerBlock);

	if (n < numBlocks * numThreadsPerBlock * m) {
		printf("The text size is too small\n");
		exit(1);
	}

	unsigned int shiftsize = determine_shiftsize(alphabet);

	//Allocate host memory for results array
	unsigned int *h_out = (unsigned int *) malloc(
			numBlocks * numThreadsPerBlock * sizeof(unsigned int));
	memset(h_out, 0, numBlocks * numThreadsPerBlock * sizeof(unsigned int));

	//Allocate 1D device memory
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_text, n * sizeof ( unsigned char ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_SHIFT, shiftsize * sizeof ( int ) ));
	checkCudaErrors(
			cudaMalloc ( ( void** ) &d_PREFIX_size, shiftsize * sizeof ( int ) ));

	//Allocate 2D device memory
	checkCudaErrors(
			cudaMallocPitch ( &d_pattern, &pattern_pitch, m * sizeof ( unsigned char ), p_size ));
	checkCudaErrors(
			cudaMallocPitch ( &d_PREFIX_value, &PREFIX_value_pitch, p_size * sizeof ( int ), shiftsize ));
	checkCudaErrors(
			cudaMallocPitch ( &d_PREFIX_index, &PREFIX_index_pitch, p_size * sizeof ( int ), shiftsize ));

	//printf("pattern_pitch %i PREFIX_value_pitch %i PREFIX_index_pitch %i\n",
	//		pattern_pitch, PREFIX_value_pitch, PREFIX_index_pitch);

	//Copy 1D host memory to device
	checkCudaErrors(
			cudaMemcpy ( d_text, text, n * sizeof ( unsigned char ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_out, h_out, numBlocks * numThreadsPerBlock * sizeof ( unsigned int ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_SHIFT, SHIFT, shiftsize * sizeof ( int ), cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy ( d_PREFIX_size, PREFIX_size, shiftsize * sizeof ( int ), cudaMemcpyHostToDevice ));

	//Copy 2D host memory to device
	checkCudaErrors(
			cudaMemcpy2D ( d_pattern, pattern_pitch, pattern, m * sizeof ( unsigned char ), m * sizeof ( unsigned char ), p_size, cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy2D ( d_PREFIX_value, PREFIX_value_pitch, PREFIX_value, p_size * sizeof ( int ), p_size * sizeof ( int ), shiftsize, cudaMemcpyHostToDevice ));
	checkCudaErrors(
			cudaMemcpy2D ( d_PREFIX_index, PREFIX_index_pitch, PREFIX_index, p_size * sizeof ( int ), p_size * sizeof ( int ), shiftsize, cudaMemcpyHostToDevice ));

	//Create timer
	cudaEvent_t start, stop;

	float time;

	//Create the timer events
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//Start the event clock	
	cudaEventRecord(start, 0);

	//Executing kernel in the device
	wm_kernel1<<<dimGrid, dimBlock>>>(d_SHIFT, d_PREFIX_value, d_PREFIX_index,
			d_PREFIX_size, d_pattern, d_text, d_out, pattern_pitch,
			PREFIX_value_pitch, PREFIX_index_pitch, m, n, p_size, alphabet,
			numBlocks);

	checkCUDAError("kernel invocation");

	cudaEventRecord(stop, 0);

	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	//Get back the results from the device
	cudaMemcpy(h_out, d_out,
			numBlocks * numThreadsPerBlock * sizeof(unsigned int),
			cudaMemcpyDeviceToHost);

	//Look at the results
	int i, matches = 0;
	for (i = 0; i < numBlocks * numThreadsPerBlock; i++)
		matches += h_out[i];

	//printf("Kernel 1 matches \t%i\t time \t%f\n", matches, time / 1000);
	*gpuTime = time / 1000;
	//Free host and device memory
	free(h_out);

	cudaFree(d_text);
	cudaFree(d_pattern);
	cudaFree(d_SHIFT);
	cudaFree(d_PREFIX_value);
	cudaFree(d_PREFIX_index);
	cudaFree(d_PREFIX_size);
	cudaFree(d_out);
	return matches;
}

