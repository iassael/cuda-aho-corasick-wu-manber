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

#include "../smatcher.h"

/// initialize the empty-table
void sbom_init ( struct sbom_table *g, int alphabet, int *state_transition ) {

	g->zerostate = NULL;
	g->patterncounter = 0;

	//Create the root note
	g->zerostate = malloc ( sizeof ( struct sbom_state ) );

	if ( !g->zerostate )
		fail ( "Could not allocate memory\n" );

	g->idcounter = 1;
	g->zerostate->id = 0;

	g->zerostate->F = NULL;
	
	//Set Supply(q_0) := fail
	g->zerostate->fail = NULL;
	
	g->zerostate->next = ( struct sbom_state ** ) malloc ( alphabet * sizeof ( struct sbom_state * ) );

	//Set all alphabet bytes of root node->next to 0
	memset ( g->zerostate->next, 0, alphabet * sizeof ( struct sbom_state * ) );
	
	//Set all cells of transition table for state 0 to 0
	int i;
	
	for ( i = 0; i < alphabet; i++ )
		state_transition[i] = 0;
}

// Insert a string to the tree
void sbom_addstring ( struct sbom_table *g, unsigned int i, unsigned char *string, int m, int p_size, int alphabet, int *state_transition, unsigned int *state_final_multi ) {
	
	struct sbom_state *state, *next = NULL, *k;
	int j, done = 0;

	// as long as next already exists follow them
	j = m - 1;
	state = g->zerostate;
	
	while ( !done && ( next = state->next[*( string + j )] ) != NULL ) {

		state = next;

		if ( j <= 0 )
			done = 1;
			
		j--;
	}

	// not done yet
	if ( !done ) {

		while ( j >= 0 ) {
			// Create new state
			next = malloc ( sizeof ( struct sbom_state ) );

			if ( !next )
				fail ( "Could not allocate memory\n" );
				
			next->next = ( struct sbom_state ** ) malloc ( alphabet * sizeof ( struct sbom_state * ) );

			next->id = g->idcounter++;
			next->F = NULL;
			
			state_transition[state->id * alphabet + *( string + j )] = next->id;
			
			//Store the pointer to the new state in an array so it can be free'ed at the end
			pointer_array[next->id - 1] = next;
			
			//printf("Created link from state %i to %i for character %i  (j = %i)\n", state->id, next->id, *( string + j ), j );

			//Set all alphabet bytes of the next node's->next to 0
			//This is the _extended_ Aho-Corasick algorithm. A complete automaton is used where all states 
			//have an outgoing transition for every alphabet character of the alphabet
			memset ( next->next, 0, alphabet * sizeof ( struct sbom_state * ) );

			state->next[*( string + j )] = next;
			
			k = state->fail;
			
			while ( k != NULL && k->next[*( string + j )] == NULL ) {
				
				k->next[*( string + j )] = next;
				
				state_transition[k->id * alphabet + *( string + j )] = next->id;
				
				//printf("	Created additional link from state %i to %i for character %i\n", k->id, next->id, *( string + j ) );
				
				k = k->fail;
			}

			if ( k != NULL )
				next->fail = k->next[*( string + j )];
			else
				next->fail = g->zerostate;
			
			state = next;
			
			j--;
		}
	}
	
	//printf("	Currently at state %i\n", state->id);
	
	//After finishing with the previous characters of the keyword, add the terminal state to F(q)
	if ( !state->F ) {
		
		//In the worst case, one state can correspond to all p_size patterns, needing p_size * number_of_terminal_states memory. A number of 200 indices should suffice.
		//state->F = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * p_size );
		state->F = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * 200 );
		
		if ( !state->F )
			fail ( "Could not allocate memory\n" );

		state->num = 0;
	}	
	
	//Add the row index to the F array
	state->F[state->num] = g->patterncounter;
	
	//printf("		Added pattern %i to F[%i] of state %i\n", g->patterncounter, state->num, state->id);
	
	//Use state_final_multi[state][0] to store the number of matching patterns, enumerate them in cells state_final_multi[state][1-200]
	state_final_multi[state->id * 200] = state->num + 1;
	state_final_multi[state->id * 200 + state->num + 1] = g->patterncounter;
	
	state->num++;
		
	g->patterncounter++;
}

unsigned int search_sbom ( unsigned char **pattern, int m, unsigned char *text, int n, struct sbom_table *table ) {

	struct sbom_state *head = table->zerostate;
	struct sbom_state *r, *s;
	
	unsigned int i;
	
	int column = m - 1, matches = 0, j;
	
	while ( column < n ) {
		
		r = head;
		j = 0;
			
		while ( j < m && ( s = r->next[*( text + column - j )] ) != NULL ) {
			
			//printf("(%i) Going from %i to %i by %i\n", column - j, r->id, s->id, *( text + column - j ));
			
			r = s;
			
			j++;
		}

		//Verify all patterns in F(q) with the input string
		if ( r->F != NULL && r->num > 0 && j == m ) {
			
			for ( i = 0; i < r->num; i++ ) {

				if ( memcmp ( pattern[r->F[i]], text + column - m + 1, m ) == 0 ) {
					matches++;
					
					//printf("match of %i %i at %i\n", r->id, r->F[i], column);

					break;
				}
			}
			
			column++;
		}
		else	
			column += MAX ( m - j, 1);
	}

	return matches;
}

struct sbom_table *preproc_sbom ( unsigned char **pattern, int m, int p_size, int alphabet, int *state_transition, unsigned int *state_final_multi ) {

	unsigned int i;
	
	struct sbom_table *table;

	// allocate memory for the table

	table = malloc ( sizeof ( struct sbom_table ) );

	if ( !table )
		fail ( "Could not initialize table\n" );

	sbom_init ( table, alphabet, state_transition );

	for ( i = 0; i < p_size; i++ )
		sbom_addstring ( table, i, pattern[i], m, p_size, alphabet, state_transition, state_final_multi );
		
	return table;
}

void free_sbom ( struct sbom_table *table, int m ) {

	int i;
	
	//We know exactly how many states we stored in the pointer_array ( table->idcounter - 1 )
	for ( i = 0; i < table->idcounter - 1; i++ )
		if ( pointer_array[i] ) {
		
			if ( pointer_array[i]->F )
				free ( pointer_array[i]->F );

			free ( pointer_array[i]->next );
			free ( pointer_array[i] );
		}

	free ( table->zerostate->next );	
	free ( table->zerostate );
	free ( table );
}

