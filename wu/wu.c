#include "../smatcher.h"

void wu_determine_shiftsize(int alphabet) {

	//the maximum size of the hash value of the B-size suffix of the patterns for the Wu-Manber algorithm
	if (alphabet == 2)
		shiftsize = 22; // 1 << 2 + 1 << 2 + 1 + 1

	else if (alphabet == 4)
		shiftsize = 64; // 3 << 2 + 3 << 2 + 3 + 1

	else if (alphabet == 8)
		shiftsize = 148; // 7 << 2 + 7 << 2 + 7 + 1

	else if (alphabet == 20)
		shiftsize = 400; // 19 << 2 + 19 << 2 + 19 + 1

	else if (alphabet == 128)
		shiftsize = 2668; // 127 << 2 + 127 << 2 + 127 + 1

	else if (alphabet == 256)
		shiftsize = 5356; //304 << 2 + 304 << 2 + 304 + 1

	else if (alphabet == 512)
		shiftsize = 10732; //560 << 2 + 560 << 2 + 560 + 1

	else if (alphabet == 1024)
		shiftsize = 21484; //1072 << 2 + 1072 << 2 + 1072 + 1

	else
		fail("The alphabet size is not supported by wu-manber\n");
}

unsigned int search_wu(unsigned char **pattern, int m, int p_size,
		unsigned char *text, int n, int *SHIFT, int *PREFIX_value,
		int *PREFIX_index, int *PREFIX_size) {

	int column = m - 1, i;

	unsigned int hash1, hash2;

	unsigned int matches = 0;

	size_t shift;

	while (column < n) {

		hash1 = text[column - 2];
		hash1 <<= m_nBitsInShift;
		hash1 += text[column - 1];
		hash1 <<= m_nBitsInShift;
		hash1 += text[column];

		shift = SHIFT[hash1];

		//printf("column %i hash1 = %i shift = %i\n", column, hash1, shift);

		if (shift == 0) {

			hash2 = text[column - m + 1];
			hash2 <<= m_nBitsInShift;
			hash2 += text[column - m + 2];

			//printf("hash2 = %i PREFIX[hash1].size = %i\n", hash2, PREFIX[hash1].size);

			//For every pattern with the same suffix as the text
			for (i = 0; i < PREFIX_size[hash1]; i++) {

				//If the prefix of the pattern matches that of the text
				if (hash2 == PREFIX_value[hash1 * p_size + i]) {

					//Compare directly the pattern with the text
					if (memcmp(pattern[PREFIX_index[hash1 * p_size + i]],
							text + column - m + 1, m) == 0) {

						matches++;

						//printf("Match of pattern index %i at %i\n", PREFIX_index[hash1 * p_size + i], column);

						break;
					}

				}
			}

			column++;
		} else
			column += shift;
	}

	return matches;
}

void preproc_wu(unsigned char **pattern, int m, int p_size, int alphabet, int B,
		int *SHIFT, int *PREFIX_value, int *PREFIX_index, int *PREFIX_size) {

	unsigned int j, q, hash;

	size_t shiftlen, prefixhash;

	for (j = 0; j < p_size; ++j) {

		//add each 3-character subpattern (similar to q-grams)
		for (q = m; q >= B; --q) {

			hash = pattern[j][q - 2 - 1]; // bring in offsets of X in pattern j
			hash <<= m_nBitsInShift;
			hash += pattern[j][q - 1 - 1];
			hash <<= m_nBitsInShift;
			hash += pattern[j][q - 1];

			//printf("hash = %i pattern[%i][%i] = %i pattern[%i][%i] = %i pattern[%i][%i] = %i\n", hash, j, q - 2 - 1, pattern[j][q - 2 - 1], j, q - 2, pattern[j][q - 2], j, q - 1, pattern[j][q - 1], j );

			shiftlen = m - q;

			SHIFT[hash] = MIN(SHIFT[hash], shiftlen);

			//calculate the hash of the prefixes for each pattern
			if (shiftlen == 0) {

				prefixhash = pattern[j][0];
				prefixhash <<= m_nBitsInShift;
				prefixhash += pattern[j][1];

				PREFIX_value[hash * p_size + PREFIX_size[hash]] = prefixhash;
				PREFIX_index[hash * p_size + PREFIX_size[hash]] = j;

				PREFIX_size[hash]++;

				//printf("%i) PREFIX[%i].value[%i] = %i PREFIX[%i].index[%i] = %i\n", j, hash, PREFIX[hash].size - 1, PREFIX[hash].value[PREFIX[hash].size - 1], hash, PREFIX[hash].size - 1, hashmap[j].index );
			}
		}
	}
}

unsigned int search_wu2(unsigned char *pattern, int m, int p_size,
		unsigned char *text, int n, int *SHIFT, int *PREFIX_value,
		int *PREFIX_index, int *PREFIX_size) {

	int column = m - 1, i;

	unsigned int hash1, hash2;

	unsigned int matches = 0;

	size_t shift;

	while (column < n) {

		hash1 = text[column - 2];
		hash1 <<= m_nBitsInShift;
		hash1 += text[column - 1];
		hash1 <<= m_nBitsInShift;
		hash1 += text[column];

		shift = SHIFT[hash1];

		//printf("column %i hash1 = %i shift = %i\n", column, hash1, shift);

		if (shift == 0) {

			hash2 = text[column - m + 1];
			hash2 <<= m_nBitsInShift;
			hash2 += text[column - m + 2];

			//printf("hash2 = %i PREFIX[hash1].size = %i\n", hash2, PREFIX[hash1].size);

			//For every pattern with the same suffix as the text
			for (i = 0; i < PREFIX_size[hash1]; i++) {

				//If the prefix of the pattern matches that of the text
				if (hash2 == PREFIX_value[hash1 * p_size + i]) {

					//Compare directly the pattern with the text
					if (memcmp(pattern + (PREFIX_index[hash1 * p_size + i] * m),
							text + column - m + 1, m) == 0) {

						matches++;

						//printf("Match of pattern index %i at %i\n", PREFIX_index[hash1 * p_size + i], column);

						break;
					}

				}
			}

			column++;
		} else
			column += shift;
	}

	return matches;
}

void preproc_wu2(unsigned char *pattern, int m, int p_size, int alphabet, int B,
		int *SHIFT, int *PREFIX_value, int *PREFIX_index, int *PREFIX_size) {

	unsigned int j, q, hash;

	size_t shiftlen, prefixhash;

	for (j = 0; j < p_size; ++j) {

		//add each 3-character subpattern (similar to q-grams)
		for (q = m; q >= B; --q) {

			hash = pattern[j * m + (q - 2 - 1)]; // bring in offsets of X in pattern j
			hash <<= m_nBitsInShift;
			hash += pattern[j * m + (q - 1 - 1)];
			hash <<= m_nBitsInShift;
			hash += pattern[j * m + (q - 1)];

			//printf("hash = %i pattern[%i][%i] = %i pattern[%i][%i] = %i pattern[%i][%i] = %i\n", hash, j, q - 2 - 1, pattern[j][q - 2 - 1], j, q - 2, pattern[j][q - 2], j, q - 1, pattern[j][q - 1], j );

			shiftlen = m - q;

			SHIFT[hash] = MIN(SHIFT[hash], shiftlen);

			//calculate the hash of the prefixes for each pattern
			if (shiftlen == 0) {

				prefixhash = pattern[j * m];
				prefixhash <<= m_nBitsInShift;
				prefixhash += pattern[j * m + 1];

				PREFIX_value[hash * p_size + PREFIX_size[hash]] = prefixhash;
				PREFIX_index[hash * p_size + PREFIX_size[hash]] = j;

				PREFIX_size[hash]++;

				//printf("%i) PREFIX[%i].value[%i] = %i PREFIX[%i].index[%i] = %i\n", j, hash, PREFIX[hash].size - 1, PREFIX[hash].value[PREFIX[hash].size - 1], hash, PREFIX[hash].size - 1, hashmap[j].index );
			}
		}
	}
}

