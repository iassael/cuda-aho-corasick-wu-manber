#include "../smatcher.h"

// A structure for holding the hash value and the pattern for the 8-byte Rabin-Karp implementation
typedef struct {

	uint32_t hs;
	uint8_t pat[32];
	int index;

} pat_hs_t32;

//Scanner that provides final matching for 8-byte patterns with Rabin-Karp.
typedef struct {
	
	// 2-level hash table
	uint8_t hs2[256*32];

	// Table holding the patterns and their hash values. This table is ordered according to the hash values
	pat_hs_t32 *patterns;
	
	// Position of the first empty slot in the pattern table
	int pos;

} sog_scanner32;

sog_scanner32 *scanner32;

#define GET32(address) (((uint32_t)((address)[0]) << 24) + ((uint32_t)((address)[1]) << 16) + ((uint32_t)((address)[2]) << 8) + (address)[3])

//Compare two patterns using their hash values
static int compSign ( const void* s1, const void* s2 ) {

	uint32_t h1 = ( (pat_hs_t32 *) s1 )->hs;
	uint32_t h2 = ( (pat_hs_t32 *) s2 )->hs;

	if (h1 < h2)
		return -1;
	else if (h1 == h2)
		return 0;
	else
		return 1;
}

int sog_rkbt_verification32 ( unsigned char *text, int m, int p_size ) {

	uint32_t hs = ( GET32((text)) ^ GET32((text + 4)) ) ^ ( GET32((text + 8)) ^ GET32((text + 12)) ) ^ ( GET32((text + 16)) ^ GET32((text + 20)) ) ^ ( GET32((text + 24)) ^ GET32((text + 28)) );
	
/*	printf("text = %c%c%c%c\n", *(text), *(text + 1), *(text + 2), *(text + 3));
	printf("text = %s\n", text);
	printf("text hs = %i\n", hs);
*/
	uint16_t hs2level = (uint16_t) ((hs >> 16) ^ hs);

	//printf("---%s\n", scanner32->patterns[lookfor].pat);
	
	/* check 2-level hash */
	if ( scanner32->hs2[hs2level >> 3] & mask[hs2level & 0x07] ) {

		int lo = 0;
		int hi = p_size - 1;
		int mid;
		uint32_t hs_pat;
			
		// do the binary search
		while ( hi >= lo ) {

			mid = ( lo + hi ) / 2;
			hs_pat = scanner32->patterns[mid].hs;
		
			//printf("	mid = %i hs = %i hs_pat = %i index = %i pat = %s \n", mid, hs, scanner32->patterns[mid].hs, scanner32->patterns[mid].index, scanner32->patterns[mid].pat);
		
			if ( hs > hs_pat )
				lo = ++mid;

			else if ( hs < hs_pat )
				hi = --mid;
			
			//if text hash equals pattern hash verify the match
			else {

				// check for duplicates and patterns with same hash
				while ( mid > 0 && hs == scanner32->patterns[mid - 1].hs )
					mid--;
				
				do {
					//printf("%c%c%c%c%c%c%c%c - %s\n", *(index - 7), *(index - 6), *(index - 5), *(index - 4), *(index - 3), *(index - 2), *(index - 1), *(index - 0), scanner32->patterns[mid].pat );
					
					if ( memcmp ( text, scanner32->patterns[mid].pat, 32 ) == 0 )
						return 1;

					mid++;
				
				} while ( mid < p_size && hs == scanner32->patterns[mid].hs );
			
				break;
			}
		}
	}
	return -1;
}

unsigned int search_sog32 ( unsigned char **pattern, int m, unsigned char *text, int n, int p_size, int B ) {

	register uint32_t E = 0xffffffff;

	int column, matches = 0;
	
	for ( column = 0; column < n - 2; column++ ) {
		
		//printf("\ncolumn %i character %c\n", column, *(text + column));
			
		//printbinary(E, 8);
		
		//printf("hs: %i T[hs]: %i ", GET3GRAM( text + column ), T8[GET3GRAM( text + column )] );
		
		//printbinary(T8[GET3GRAM( text + column )], 8);
	
		E = (E << 1) | T8[GET3GRAM( text + column )];
		
		//printbinary(E, 8);
		
		//printbinary(E & 0x20, 8);
		
		if ( E & 0x20000000 )
			continue;
			
		//printf("potential match at column %i\n", column + B - 1);
		
		if ( sog_rkbt_verification32 ( (unsigned char *)text + column - m + B, m, p_size ) != -1 )
			matches++;
		
	}

	return matches;
}

static void sog_add_pattern2 ( uint8_t *pattern, int m, int p_size ) {

	int i;

	uint32_t hs;
	uint16_t hs2level;

	if ( scanner32->pos < p_size ) {

		//add pattern
		for ( i = 0; i < m; i++ )
			scanner32->patterns[scanner32->pos].pat[i] = pattern[i];
			
		//add index
		scanner32->patterns[scanner32->pos].index = scanner32->pos;

		// Count hash
		scanner32->patterns[scanner32->pos].hs = ( GET32(pattern) ^ GET32(&pattern[4]) ) ^ ( GET32(&pattern[8]) ^ GET32(&pattern[12]) ) ^ ( GET32(&pattern[16]) ^ GET32(&pattern[20]) ) ^ ( GET32(&pattern[24]) ^ GET32(&pattern[28]) );
		
		//printf("scanner32->patterns[%i].hs = %i\n", scanner32->pos, scanner32->patterns[scanner32->pos].hs);

		// Count 2-level hash
		hs = scanner32->patterns[scanner32->pos].hs;
		hs2level = ( uint16_t ) ( ( hs >> 16 ) ^ hs );
		
		scanner32->hs2[hs2level >> 3] |= mask[hs2level & 0x07];
		
		scanner32->pos++;
	}
}

static void sog_add_pattern ( uint8_t *pattern, int m, int p_size ) {
	
	uint8_t *index = &pattern[0];
	uint8_t *limit = &pattern[31];
	
	unsigned int i = 0;

	uint32_t hs;
	
	sog_add_pattern2 ( pattern, m, p_size );

	while ( index < limit ) {
		hs = GET3GRAM( index );
		
		//printbinary(hs, 32);
		//printf("hs: %i T[hs]: %i ", hs, T[hs]);
		
		T32[hs] &= 0xffffffff - ( 1 << i );

		//printbinary(T[hs], 8);
		
		index++;
		i++;
	}
	
	//printf("\n");
}

static void sog_reset_patterns ( int m ) {
	
	unsigned int i;

	for ( i = 0; i < SIZE_3GRAM_TABLE; i++ )
		T32[i] = 0xffffffff;

	scanner32->pos = 0;

	// Reset 2-level hashes
	for ( i = 0; i < 32 * 256; i++ )
		scanner32->hs2[i] = 0x00;
}

void sog_init32 ( int p_size) {

	scanner32 = malloc ( sizeof ( sog_scanner32 ) );
	scanner32->patterns = malloc ( p_size * sizeof ( pat_hs_t32 ) );
}

void sog_free32 () {

	free ( scanner32->patterns );
	free ( scanner32 );
}

void preproc_sog32 ( unsigned char **pattern, int m, int p_size ) {

	unsigned int i;
	
	sog_reset_patterns ( p_size );
	
	for ( i = 0; i < p_size; i++ )
		sog_add_pattern ( pattern[i], m, p_size );

	//Sort the patterns so that binary search can be used
	qsort ( scanner32->patterns, p_size, sizeof( pat_hs_t32 ), compSign );
}

