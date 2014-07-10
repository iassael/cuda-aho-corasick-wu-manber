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

#include <mpi.h>

#include "smatcher.h"

void usage(void) {

	printf(
			"Smatcher - Sequential two dimensional pattern matching toolkit using multiple pattern matching algorithms.\n");
	printf(
			"Usage: <multiple algorithm> <2d algorithm> -m <m> -p_size <p_size> -n <n> -alphabet <alphabet>\n");
	printf("-h,--help,\t\t print this help message\n");
	printf("-c,\t\t create the data files\n");

	exit(0);
}

void select_data_file(unsigned int m, unsigned int n, unsigned int alphabet,
		char *pattern_filename, char *text_filename, int create_data) {

	sprintf(pattern_filename, "../data-cuda-multi/pattern/%i/%i/%i/pattern", n,
			m, alphabet);

	switch (n) {
	case 3999744:
		if ((alphabet != 2 && alphabet != 8))
			fail("For random texts, you must use an alphabet size of 2 or 8\n");

		if (alphabet == 2)
			sprintf(text_filename, "../data-cuda-multi/text/text2");
		else if (alphabet == 8)
			sprintf(text_filename, "../data-cuda-multi/text/text8");

		if (create_data)
			create_multiple_pattern_with_hits(m, n, 100000, text_filename,
					pattern_filename);

		break;

	case 1903104:
		if (alphabet != 128)
			fail("For english text, you must use an alphabet size of 128\n");

		sprintf(text_filename, "../data-cuda-multi/text/world192.txt");

		if (create_data)
			create_multiple_pattern_with_hits(m, n, 100000, text_filename,
					pattern_filename);

		break;
		//287x16128
	case 4628736:
		if (alphabet != 4)
			fail("For DNA sequences, you must use an alphabet size of 4\n");

		sprintf(text_filename, "../data-cuda-multi/text/E.coli2");

		if (create_data)
			create_multiple_pattern_with_hits(m, n, 100000, text_filename,
					pattern_filename);

		break;
	case 177649920:
		if (alphabet != 20)
			fail("For swiss-prot, you must use an alphabet size of 20\n");

		sprintf(text_filename, "../data-cuda-multi/text/swiss-prot");

		if (create_data)
			create_multiple_pattern_with_hits(m, n, 100000, text_filename,
					pattern_filename);

		break;
	case 10821888:
		if (alphabet != 20)
			fail("For A_thaliana.faa, you must use an alphabet size of 20\n");

		sprintf(text_filename, "../data-cuda-multi/text/A_thaliana.faa");

		if (create_data)
			create_multiple_pattern_with_hits(m, n, 100000, text_filename,
					pattern_filename);

		break;
	case 116234496:
		if (alphabet != 4)
			fail("For A_thaliana.fna, you must use an alphabet size of 4\n");

		sprintf(text_filename, "../data-cuda-multi/text/A_thaliana.fna");

		if (create_data)
			create_multiple_pattern_with_hits(m, n, 100000, text_filename,
					pattern_filename);

		break;

	case 100:
		if (alphabet != 2)
			fail("The debug text uses a binary alphabet\n");

		sprintf(text_filename, "../data-cuda-multi/text/debug");
		sprintf(pattern_filename, "../data-cuda-multi/pattern/debug");

		break;
	default:
		fail("Please select an appropriate text size\n");
		break;
	}
}

/*void multiac(unsigned char **pattern, int m, unsigned char *text, int n,
 int p_size, int alphabet, unsigned int *state_transition,
 unsigned int *state_supply, unsigned int *state_final) {

 double t1, t2, t3, preproc_time = 0, search_time = 0, running_time = 0;
 int matches;

 struct ac_table *table;

 //for ( i = 0; i < 10; i++ ) {

 t1 = MPI_Wtime();

 table = preproc_ac(pattern, m, p_size, alphabet, state_transition,
 state_supply, state_final);

 t2 = MPI_Wtime();

 matches = search_ac(text, n, table);

 t3 = MPI_Wtime();

 preproc_time += (t2 - t1);

 search_time += (t3 - t2);

 running_time += (t3 - t1);

 free_ac(table, alphabet);
 //}

 printf("search_ac matches \t%i\t time \t%f\n", matches, search_time);
 }

 void multish(unsigned char **pattern, int m, unsigned char *text, int n,
 int p_size, int alphabet, unsigned int *state_transition,
 unsigned int *state_final, int *bmBc) {

 double t1, t2, t3, preproc_time = 0, search_time = 0, running_time = 0;
 int matches;

 struct ac_table *table;

 //for ( i = 0; i < 10; i++ ) {

 t1 = MPI_Wtime();

 //Preprocessing using the Horspool algorithm
 preBmBc(pattern, m, p_size, alphabet, bmBc);

 //Creating the Set Horspool automaton AND translating the pattern into an 1D array of symbols
 table = preproc_sh(pattern, m, p_size, alphabet, state_transition,
 state_final);

 t2 = MPI_Wtime();

 matches = search_sh(m, text, n, table, bmBc);

 t3 = MPI_Wtime();

 preproc_time += (t2 - t1);

 search_time += (t3 - t2);

 running_time += (t3 - t1);

 free_sh(table, alphabet);
 //}

 printf("search_sh matches \t%i\t time \t%f\n", matches, search_time);
 }

 void multisbom(unsigned char **pattern, int m, unsigned char *text, int n,
 int p_size, int alphabet, unsigned int *state_transition,
 unsigned int *state_final_multi) {

 double t1, t2, t3, preproc_time = 0, search_time = 0, running_time = 0;
 int matches;

 struct sbom_table *table;

 //for ( i = 0; i < 10; i++ ) {

 pointer_array = malloc(p_size * m * sizeof(struct sbom_state));

 t1 = MPI_Wtime();

 table = preproc_sbom(pattern, m, p_size, alphabet, state_transition,
 state_final_multi);

 t2 = MPI_Wtime();

 matches = search_sbom(pattern, m, text, n, table);

 t3 = MPI_Wtime();

 preproc_time += (t2 - t1);

 search_time += (t3 - t2);

 running_time += (t3 - t1);

 free_sbom(table, m);

 free(pointer_array);
 //}

 printf("search_sbom matches \t%i\t time \t%f\n", matches, search_time);
 }

 //void multiwm ( unsigned char **pattern, int m, unsigned char *text, int n, int p_size, int alphabet, int B, int *SHIFT, struct prefixArray *PREFIX ) {
 */
void multiwm(unsigned char **pattern, int m, unsigned char *text, int n,
		int p_size, int alphabet, int B, int *SHIFT, int *PREFIX_value,
		int *PREFIX_index, int *PREFIX_size) {

	double t1, t2, t3, preproc_time = 0, search_time = 0, running_time = 0;
	int matches;

	//for ( i = 0; i < 10; i++ ) {

	t1 = MPI_Wtime();

	preproc_wu(pattern, m, p_size, alphabet, B, SHIFT, PREFIX_value,
			PREFIX_index, PREFIX_size);

	t2 = MPI_Wtime();

	matches = search_wu(pattern, m, p_size, text, n, SHIFT, PREFIX_value,
			PREFIX_index, PREFIX_size);

	t3 = MPI_Wtime();

	preproc_time += (t2 - t1);

	search_time += (t3 - t2);

	running_time += (t3 - t1);

	//}

	printf("search_wm matches \t%i\t time \t%f\n", matches, search_time);
}
void multiwm2(unsigned char *pattern, int m, unsigned char *text, int n,
		int p_size, int alphabet, int B, int *SHIFT, int *PREFIX_value,
		int *PREFIX_index, int *PREFIX_size) {

	double t1, t2, t3, preproc_time = 0, search_time = 0, running_time = 0;
	int matches;

	//for ( i = 0; i < 10; i++ ) {

	t1 = MPI_Wtime();

	preproc_wu2(pattern, m, p_size, alphabet, B, SHIFT, PREFIX_value,
			PREFIX_index, PREFIX_size);

	t2 = MPI_Wtime();

	matches = search_wu2(pattern, m, p_size, text, n, SHIFT, PREFIX_value,
			PREFIX_index, PREFIX_size);

	t3 = MPI_Wtime();

	preproc_time += (t2 - t1);

	search_time += (t3 - t2);

	running_time += (t3 - t1);

	//}

	printf("search_wm2 matches \t%i\t time \t%f\n", matches, search_time);
}

/*void multisog(uint8_t **T8, uint32_t **scanner_hs, int **scanner_index,
 uint8_t **scanner_hs2, unsigned char **pattern, int m,
 unsigned char *text, int n, int p_size, int alphabet, int B) {

 double t1, t2, t3, preproc_time = 0, search_time = 0, running_time = 0;
 int i, matches;

 t1 = MPI_Wtime();
 preproc_sog8(T8, scanner_hs, scanner_index, scanner_hs2, pattern, m, text,
 n, p_size, B);
 t2 = MPI_Wtime();
 matches = search_sog8(T8, scanner_hs, scanner_index, scanner_hs2, pattern,
 m, text, n, p_size, B);
 t3 = MPI_Wtime();

 preproc_time += (t2 - t1);

 search_time += (t3 - t2);

 running_time += (t3 - t1);

 printf("search_sog matches \t%i\t time \t%f\n", matches, search_time);
 }*/

int main(int argc, char **argv) {

	// Initialize MPI state
	MPI_Init(&argc, &argv);

	// Get our MPI node number and node count
	int commSize, commRank;
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
	MPI_Barrier(MPI_COMM_WORLD);

	int m = 0, p_size = 0, n = 0, nFull = 0, alphabet = 0, B = 3, create_data =
			0;

	double timeScatter = 0, timeGather = 0, timeReadFile = 0,
			timeExecuteCPU = 0;

	int i, j;

	//char text_filename[100], pattern_filename[100];
	char *text_filename = (char *) malloc(100 * sizeof(char));
	char *pattern_filename = (char *) malloc(100 * sizeof(char));

	//Scan command line arguments
	for (i = 1; i < argc; i++) {

		if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0)
			usage();

		if (strcmp(argv[i], "-m") == 0)
			m = atoi(argv[i + 1]);

		if (strcmp(argv[i], "-n") == 0)
			nFull = atoi(argv[i + 1]);

		if (strcmp(argv[i], "-p_size") == 0)
			p_size = atoi(argv[i + 1]);

		if (strcmp(argv[i], "-alphabet") == 0)
			alphabet = atoi(argv[i + 1]);

		if (strcmp(argv[i], "-c") == 0)
			create_data = 1;
	}

	if (m == 0 || nFull == 0 || p_size == 0 || alphabet == 0)
		usage();

	if (p_size > 100000)
		fail("Only up to 100.000 patterns are supported\n");

	if (commSize > 1)
		n = (float) ceil((double) nFull / (float) commSize) + (m - 1);
	else
		n = nFull;

	if (commRank == 0) {
		//Determine path of pattern and text based on n
		select_data_file(m, nFull, alphabet, pattern_filename, text_filename,
				create_data);
	}

	unsigned char *textFull = (unsigned char *) malloc(
			nFull * sizeof(unsigned char));
	unsigned char *text = (unsigned char *) malloc(n * sizeof(unsigned char));

	if (text == NULL )
		fail("Failed to allocate array\n");

	unsigned char **pattern = (unsigned char **) malloc(
			p_size * sizeof(unsigned char *));

	if (pattern == NULL )
		fail("Failed to allocate array!\n");

	for (i = 0; i < p_size; i++) {
		pattern[i] = (unsigned char *) malloc(m * sizeof(unsigned char));

		if (pattern[i] == NULL )
			fail("Failed to allocate array!\n");
	}

	unsigned char *pattern2 = (unsigned char *) malloc(
			m * p_size * sizeof(unsigned char));

	//There can be a maximum of m * p_size + 1 (for the initial state) states
	int *state_transition = (int *) malloc(
			(m * p_size + 1) * alphabet * sizeof(int));
	memset(state_transition, -1, (m * p_size + 1) * alphabet * sizeof(int));

	unsigned int *state_supply = (unsigned int *) malloc(
			(m * p_size + 1) * sizeof(unsigned int));
	memset(state_supply, 0, (m * p_size + 1) * sizeof(unsigned int));

	unsigned int *state_final = (unsigned int *) malloc(
			(m * p_size + 1) * sizeof(unsigned int));
	memset(state_final, 0, (m * p_size + 1) * sizeof(unsigned int));

	//state_final_multi is reuqired for final states that correspond to more than one pattern
	unsigned int *state_final_multi = (unsigned int *) malloc(
			(m * p_size + 1) * 200 * sizeof(unsigned int));
	memset(state_final_multi, 0, (m * p_size + 1) * 200 * sizeof(unsigned int));

	int *bmBc = (int *) malloc(alphabet * sizeof(int));

	wu_determine_shiftsize(alphabet);

	m_nBitsInShift = 2;

	int *SHIFT = (int *) malloc(shiftsize * sizeof(int));

	//The hash value of the B'-character prefix of a pattern
	int *PREFIX_value = (int *) malloc(shiftsize * p_size * sizeof(int));

	//The pattern number
	int *PREFIX_index = (int *) malloc(shiftsize * p_size * sizeof(int));

	//How many patterns with the same prefix hash exist
	int *PREFIX_size = (int *) malloc(shiftsize * sizeof(int));

	for (i = 0; i < shiftsize; i++) {

		//*( *SHIFT + i ) = m - B + 1;
		SHIFT[i] = m - B + 1;
		PREFIX_size[i] = 0;
	}

	if (commRank == 0) {
		timeReadFile = MPI_Wtime();
		load_files(pattern, textFull, m, nFull, pattern_filename, text_filename,
				p_size);

		//Copy the 2d pattern array to 1d pattern2
		for (j = 0; j < p_size; j++)
			for (i = 0; i < m; i++)
				pattern2[j * m + i] = pattern[j][i];
		timeReadFile = MPI_Wtime() - timeReadFile;
	}

	//Calculate Displacement
	int* displs = (int *) malloc(commSize * sizeof(int));
	int* scounts = (int *) malloc(commSize * sizeof(int));

	for (i = 0; i < commSize; ++i) {

		int start = i * (float) ceil((double) nFull / (float) commSize);
		int stop = (i + 1) * (float) ceil((double) nFull / (float) commSize)
				+ (m - 1);
		if (stop > nFull)
			stop = nFull;

		displs[i] = start;
		scounts[i] = stop - start;
	}

// print calculated send counts and displacements for each process
//	if (0 == commRank)
//		for (i = 0; i < commSize; i++)
//			printf("scounts[%d]\t%d\tdispls[%d]\t%d\n", i, scounts[i], i,
//					displs[i]);

//Scatterv

	timeScatter = MPI_Wtime();
	MPI_Scatterv(textFull, scounts, displs, MPI_CHAR, text, n, MPI_CHAR, 0,
			MPI_COMM_WORLD);

//Broadcast pattern
	MPI_Bcast(pattern2, p_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	timeScatter = MPI_Wtime() - timeScatter;

	/*uint8_t *T8 = (uint8_t *) malloc(SIZE_3GRAM_TABLE * sizeof(uint8_t));

	 if (T8 == NULL )
	 fail("Failed to allocate array\n");

	 //holds the hash value for the 8-byte Rabin-Karp implementation
	 uint32_t *scanner_hs = (uint32_t *) malloc(p_size * sizeof(uint32_t));

	 if (scanner_hs == NULL )
	 fail("Failed to allocate array\n");

	 //holds the pattern index for the 8-byte Rabin-Karp implementation
	 int *scanner_index = (int *) malloc(p_size * sizeof(int));

	 if (scanner_index == NULL )
	 fail("Failed to allocate array\n");

	 // 2-level hash table
	 uint8_t *scanner_hs2 = (uint8_t *) malloc(256 * 32 * sizeof(uint8_t));

	 if (scanner_hs2 == NULL )
	 fail("Failed to allocate array\n");*/

//////////////////////////Serial execution//////////////////////////
	/*if (strcmp(argv[1], "ac") == 0)
	 multiac(pattern, m, text, n, p_size, alphabet, state_transition,
	 state_supply, state_final);

	 if (strcmp(argv[1], "sh") == 0)
	 multish(pattern, m, text, n, p_size, alphabet, state_transition,
	 state_final, bmBc);

	 if (strcmp(argv[1], "sbom") == 0)
	 multisbom(pattern, m, text, n, p_size, alphabet, state_transition,
	 state_final_multi);

	 if (strcmp(argv[1], "wm") == 0)*/
	MPI_Barrier(MPI_COMM_WORLD);
	timeExecuteCPU = MPI_Wtime();
	multiwm2(pattern2, m, text, n, p_size, alphabet, B, SHIFT, PREFIX_value,
			PREFIX_index, PREFIX_size);
	timeExecuteCPU = MPI_Wtime() - timeExecuteCPU;

//multiwm(pattern, m, text, n, p_size, alphabet, B, SHIFT, PREFIX_value,
//		PREFIX_index, PREFIX_size);

	/*if (strcmp(argv[1], "sog") == 0) {

	 if (m != 8)
	 fail("Set the pattern size for SOG to m = 8\n");

	 multisog(T8, scanner_hs, scanner_index, scanner_hs2, pattern, m, text,
	 n, p_size, alphabet, B);
	 }*/
////////////////////////Parallel execution////////////////////////
	/*for ( j = 0; j < m  * p_size + 1; j++ ) {
	 for ( i = 0; i < alphabet; i++ )
	 printf("%i\t", state_transition[j*alphabet + i]);
	 printf("\n");
	 }

	 for ( i = 0; i < ( m  * p_size + 1 ); i++ )
	 printf("supply[%i] = %i\n", i, state_supply[i]);

	 for ( i = 0; i < ( m  * p_size + 1 ); i++ )
	 if ( state_final[i] )
	 printf("final[%i] = %i\n", i, state_final[i]);

	 for ( j = 0; j < m  * p_size + 1; j++ ) {

	 if ( state_final_multi[j * 200] < 1 )
	 continue;

	 printf("%i: ", j);



	 for ( i = 0; i < state_final_multi[j * 200] + 1; i++ ) {
	 printf("%i ", state_final_multi[j * 200 + i]);

	 if ( state_final_multi[j * 200 + i] == 0 )
	 break;
	 }

	 printf("\n");
	 }*/

	/*if (strcmp(argv[1], "ac") == 0) {
	 cuda_ac1(m, text, n, p_size, alphabet, state_transition, state_supply,
	 state_final);
	 cuda_ac2(m, text, n, p_size, alphabet, state_transition, state_supply,
	 state_final);
	 cuda_ac3(m, text, n, p_size, alphabet, state_transition, state_supply,
	 state_final);
	 cuda_ac4(m, text, n, p_size, alphabet, state_transition, state_supply,
	 state_final);
	 cuda_ac5(m, text, n, p_size, alphabet, state_transition, state_supply,
	 state_final);
	 }

	 if (strcmp(argv[1], "sh") == 0) {
	 cuda_sh1(m, text, n, p_size, alphabet, state_transition, state_final,
	 bmBc);
	 cuda_sh2(m, text, n, p_size, alphabet, state_transition, state_final,
	 bmBc);
	 cuda_sh3(m, text, n, p_size, alphabet, state_transition, state_final,
	 bmBc);
	 cuda_sh4(m, text, n, p_size, alphabet, state_transition, state_final,
	 bmBc);
	 cuda_sh5(m, text, n, p_size, alphabet, state_transition, state_final,
	 bmBc);
	 }

	 if (strcmp(argv[1], "sbom") == 0) {
	 cuda_sbom1(pattern2, m, text, n, p_size, alphabet, state_transition,
	 state_final_multi);
	 cuda_sbom2(pattern2, m, text, n, p_size, alphabet, state_transition,
	 state_final_multi);
	 cuda_sbom3(pattern2, m, text, n, p_size, alphabet, state_transition,
	 state_final_multi);
	 cuda_sbom4(pattern2, m, text, n, p_size, alphabet, state_transition,
	 state_final_multi);
	 cuda_sbom5(pattern2, m, text, n, p_size, alphabet, state_transition,
	 state_final_multi);
	 }

	 if (strcmp(argv[1], "wm") == 0) {
	 */
	int results[5] = { 0 };
	double gpuTime[5] = { 0.0 };
	double gpuTimeMPI[5] = { 0.0 };
	double gpuTime_sum[5] = { 0.0 };
	int result = 0;

	gpuTimeMPI[0] = MPI_Wtime();
	results[0] = cuda_wm1(pattern2, m, text, n, p_size, alphabet, B, SHIFT,
			PREFIX_value, PREFIX_index, PREFIX_size, &gpuTime[0]);
	gpuTimeMPI[0] = MPI_Wtime() - gpuTimeMPI[0];
	gpuTimeMPI[1] = MPI_Wtime();
	results[1] = cuda_wm2(pattern2, m, text, n, p_size, alphabet, B, SHIFT,
			PREFIX_value, PREFIX_index, PREFIX_size, &gpuTime[1]);
	gpuTimeMPI[1] = MPI_Wtime() - gpuTimeMPI[1];
	gpuTimeMPI[2] = MPI_Wtime();
	results[2] = cuda_wm3(pattern2, m, text, n, p_size, alphabet, B, SHIFT,
			PREFIX_value, PREFIX_index, PREFIX_size, &gpuTime[2]);
	gpuTimeMPI[2] = MPI_Wtime() - gpuTimeMPI[2];
	gpuTimeMPI[3] = MPI_Wtime();
	results[3] = cuda_wm4(pattern2, m, text, n, p_size, alphabet, B, SHIFT,
			PREFIX_value, PREFIX_index, PREFIX_size, &gpuTime[3]);
	gpuTimeMPI[3] = MPI_Wtime() - gpuTimeMPI[3];
	gpuTimeMPI[4] = MPI_Wtime();
	results[4] = cuda_wm5(pattern2, m, text, n, p_size, alphabet, B, SHIFT,
			PREFIX_value, PREFIX_index, PREFIX_size, &gpuTime[4]);
	gpuTimeMPI[4] = MPI_Wtime() - gpuTimeMPI[4];

	result = results[0];
	int result_sum = 0;

// compute global sum
	MPI_Barrier(MPI_COMM_WORLD);
	timeGather = MPI_Wtime();
	MPI_Reduce(&result, &result_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	timeGather = MPI_Wtime() - timeGather;

	for (i = 0; i < 5; i++)
		MPI_Reduce(&gpuTime[i], &gpuTime_sum[i], 1, MPI_DOUBLE, MPI_SUM, 0,
				MPI_COMM_WORLD);

	if (0 == commRank) {
		printf("Total results: %d.\n", result_sum);
		printf("timeReadFile: %f.\n", timeReadFile);
		printf("timeScatter: %f.\n", timeScatter);
		printf("timeExecuteCPU: %f.\n", timeExecuteCPU);
		printf("timeGather: %f.\n", timeGather);
		for (i = 0; i < 5; i++)
			printf("gpuTime[%d]: %f.\n", (i + 1), gpuTime_sum[i] / commSize);
	}
	/*}

	 if (strcmp(argv[1], "sog") == 0) {
	 cuda_sog1(T8, scanner_hs, scanner_index, scanner_hs2, pattern2, m, text,
	 n, p_size, B);
	 cuda_sog2(T8, scanner_hs, scanner_index, scanner_hs2, pattern2, m, text,
	 n, p_size, B);
	 cuda_sog3(T8, scanner_hs, scanner_index, scanner_hs2, pattern2, m, text,
	 n, p_size, B);
	 cuda_sog4(T8, scanner_hs, scanner_index, scanner_hs2, pattern2, m, text,
	 n, p_size, B);
	 cuda_sog5(T8, scanner_hs, scanner_index, scanner_hs2, pattern2, m, text,
	 n, p_size, B);
	 }*/

//////////////////////////clean up phase//////////////////////////
// finalize MPI
	MPI_Finalize();

	free(text);
	free(textFull);

	for (i = 0; i < p_size; i++)
		free(pattern[i]);

	free(pattern);
	free(pattern2);

	free(state_transition);
	free(state_supply);
	free(state_final);
	free(state_final_multi);
	free(bmBc);

	free(SHIFT);

	free(PREFIX_value);
	free(PREFIX_index);
	free(PREFIX_size);

	/*free(T8);
	 free(scanner_hs);
	 free(scanner_index);
	 free(scanner_hs2);*/

	return 0;
}

