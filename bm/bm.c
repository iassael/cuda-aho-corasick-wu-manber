#include "../smatcher.h"

//Bad character shift
void preBmBc ( unsigned char **pattern, int m, int p_size, int alphabet, int *bmBc ) {

	unsigned int i, j;
	
	for (i = 0; i < alphabet; ++i)
		bmBc[i] = m;

	for ( j = 0; j < p_size; j++ )
		for (i = 0; i < m - 1; ++i)
			bmBc[pattern[j][i]] = MIN ( m - i - 1, bmBc[pattern[j][i]]);
}
/*
void suffixes ( unsigned char *x, int m, int *suff ) {

	int f, g, i;

	suff[m - 1] = m;
	
	printf("suff[%i] = %i\n", m - 1, suff[m - 1]);

	g = m - 1;

	for ( i = m - 2; i >= 0; --i ) {
	
		//printf("i = %i |>| g = %i AND suff[%i] = %i |<| %i\n", i, g, i + m - 1 - f, suff[i + m - 1 - f], i - g);
		
		if ( i > g && suff[i + m - 1 - f] < i - g ) {
			suff[i] = suff[i + m - 1 - f];
			
			printf("suff[%i] = suff[%i] = %i\n", i, i + m - 1 - f, suff[i]);
		}

		else {
			if ( i < g )
				g = i;
			
			f = i;
			
			while (g >= 0 && x[g] == x[g + m - 1 - f])
				--g;
			
			suff[i] = f - g;
			
			printf("suff[%i] = %i\n", i, suff[i]);
		}
	}
}

//Good suffix shift
void preBmGs( unsigned char **pattern, int m, int bmGs[] ) {

	int i, j, suff[m];

	//suffixes( pattern, m, suff );
	
	suffixes("AACAA", m, suff );
	
	printf("\n");

	for ( i = 0; i < m; ++i )
		bmGs[i] = m;

	j = 0;

	for ( i = m - 1; i >= 0; --i )
		if ( suff[i] == i + 1 )
			for ( ; j < m - 1 - i; ++j )
				if ( bmGs[j] == m )
					bmGs[j] = m - 1 - i;
					
	for ( i = 0; i < m; i++ )
		printf("bmGs[%i] = %i\n", i, bmGs[i]);
	printf("\n");

	for (i = 0; i <= m - 2; ++i)
		bmGs[m - 1 - suff[i]] = m - 1 - i;
		
	for ( i = 0; i < m; i++ )
		printf("bmGs[%i] = %i\n", i, bmGs[i]);
		
	exit(0);
}
*/
