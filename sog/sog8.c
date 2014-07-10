#include "../smatcher.h"

#define GET32(address) (((uint32_t)((address)[0]) << 24) + ((uint32_t)((address)[1]) << 16) + ((uint32_t)((address)[2]) << 8) + (address)[3])

//quicksort implementation
void swap ( int *a, int *b ) {

	int t = *a;
	*a = *b;
	*b = t;
}

void my_sort ( uint32_t *hs, int *index, int beg, int end ) {

	if ( end > beg + 1 ) {
	
		int piv = hs[beg], l = beg + 1, r = end;

		while ( l < r ) {

			if ( hs[l] <= piv )
				l++;
			else {
        			swap ( &hs[l], &hs[--r]);
	        		swap ( &index[l], &index[r]);
        		}
		}
		
		swap ( &hs[--l], &hs[beg]);
		swap ( &index[l], &index[beg]);
		my_sort ( hs, index, beg, l );
		my_sort ( hs, index, r, end );
	}
}

int sog_rkbt_verification8 ( uint32_t *scanner_hs, int *scanner_index, uint8_t *scanner_hs2, unsigned char **pattern, unsigned char *text, int m, int p_size ) {

	uint32_t hs = GET32((text)) ^ GET32((text + 4));
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
					if ( memcmp ( text, pattern[scanner_index[mid]], 8 ) == 0 )
						return 1;

					mid++;
				
				} while ( mid < p_size && hs == scanner_hs[mid] );
			
				break;
			}
		}
	}
	return -1;
}

unsigned int search_sog8 ( uint8_t *T8, uint32_t *scanner_hs, int *scanner_index, uint8_t *scanner_hs2, unsigned char **pattern, int m, unsigned char *text, int n, int p_size, int B ) {

	register uint8_t E = 0xff;

	int column, matches = 0;
	
	for ( column = 0; column < n - 2; column++ ) {
		
		E = (E << 1) | T8[GET3GRAM( text + column )];
		
		if ( E & 0x20 )
			continue;

		if ( sog_rkbt_verification8 ( scanner_hs, scanner_index, scanner_hs2, pattern, (unsigned char *)text + column - m + B, m, p_size ) != -1 )
			matches++;
	}

	return matches;
}

static void sog_add_pattern ( uint8_t **T8, int *scanner_pos, uint32_t **scanner_hs, int **scanner_index, uint8_t **scanner_hs2, uint8_t *pattern, int m, int p_size ) {
	
	uint8_t *index = &pattern[0];
	uint8_t *limit = &pattern[6];
	
	unsigned int i = 0;

	uint32_t hs, hs2;
	uint16_t hs2level;
			
	//add index
	*( *scanner_index + *scanner_pos ) = *scanner_pos;

	// Count hash
	*( *scanner_hs + *scanner_pos ) = GET32(pattern) ^ GET32(&pattern[4]);

	// Count 2-level hash
	hs2 = *( *scanner_hs + *scanner_pos );
	hs2level = ( uint16_t ) ( ( hs >> 16 ) ^ hs2 );
		
	*( *scanner_hs2 + ( hs2level >> 3 ) ) |= mask[hs2level & 0x07];
	*scanner_pos = *scanner_pos + 1;

	while ( index < limit ) {
		hs = GET3GRAM( index );
		
		*( *T8 + hs) &= 0xff - ( 1 << i );
		
		index++;
		i++;
	}
}

static void sog_reset_patterns ( uint8_t **T8, uint8_t **scanner_hs2) {
	
	unsigned int i;

	for ( i = 0; i < SIZE_3GRAM_TABLE; i++ )
		*( *T8 + i ) = 0xff;

	// Reset 2-level hashes
	for ( i = 0; i < 32 * 256; i++ )
		*( *scanner_hs2 + i ) = 0x00;
}

void preproc_sog8 ( uint8_t *T8, uint32_t *scanner_hs, int *scanner_index, uint8_t *scanner_hs2, unsigned char **pattern, int m, unsigned char *text, int n, int p_size, int B ) {

	int i;
	
	int scanner_pos = 0;

	sog_reset_patterns ( &T8, &scanner_hs2 );
		
	for ( i = 0; i < p_size; i++ )
		sog_add_pattern ( &T8, &scanner_pos, &scanner_hs, &scanner_index, &scanner_hs2, pattern[i], m, p_size );

	my_sort ( scanner_hs, scanner_index, 0, p_size );
}

