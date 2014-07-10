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

#include "list.h"

/// free an AC table from a given startnode (recursively)
void ac_free ( struct ac_state *state, int alphabet ) {

	int i;

	for ( i = 0; i < alphabet; i++ )
		if ( state->next[i] )
			ac_free ( state->next[i], alphabet );

	if ( state->output )
		free ( state->output );
		
	free ( state->next );
	free ( state );
}

/// initialize the empty-table
void ac_init ( struct ac_table *g, int alphabet, int *state_transition ) {

	g->zerostate = NULL;
	g->patterncounter = 0;

	//Create the root note
	g->zerostate = malloc ( sizeof ( struct ac_state ) );

	if ( !g->zerostate )
		fail ( "Could not allocate memory\n" );

	g->idcounter = 1;
	g->zerostate->id = 0;

	g->zerostate->output = NULL;
	
	g->zerostate->next = ( struct ac_state ** ) malloc ( alphabet * sizeof ( struct ac_state * ) );

	//Set all alphabet bytes of root node->next to 0
	memset ( g->zerostate->next, 0, alphabet * sizeof ( struct ac_state * ) );
	
	//Set all cells of transition table for state 0 to 0
	int i;
	
	for ( i = 0; i < alphabet; i++ )
		state_transition[i] = 0;
}

/// free an entire AC table
void ac_destroy ( struct ac_table *in, int alphabet ) {

	int i;

	for ( i = 0; i < alphabet; i++ )
		if ( in->zerostate->next[i] && in->zerostate->next[i]->id > 0 ) {
			ac_free ( in->zerostate->next[i], alphabet );
			in->zerostate->next[i] = NULL;
		}
	free ( in->zerostate->next );
	free ( in->zerostate );
}

void ac_maketree ( struct ac_table *g, int alphabet, unsigned int *state_supply ) {

	struct list *list = NULL;
	struct ac_state *state, *s, *cur;
	int i/*, j*/;

	// Set all NULL transitions of 0 state to point to itself
	for ( i = 0; i < alphabet; i++ ) {
		if ( !g->zerostate->next[i] )
			g->zerostate->next[i] = g->zerostate;
		else {
			list = list_append ( list, g->zerostate->next[i] );
			g->zerostate->next[i]->fail = g->zerostate;
		}
	}

	// Set fail() for depth > 0
	while ( list ) {
		
		cur = ( struct ac_state * )list->id;
		
		for ( i = 0; i < alphabet; i++ ) {
			
			s = cur->next[i];
			
			if ( s ) {
				
				list = list_append ( list, s );
				state = cur->fail;
				
				while ( !state->next[i] )
					state = state->fail;
				
				s->fail = state->next[i];
				
				state_supply[s->id] = s->fail->id;

				//printf("Created additional link from state %i to state %i\n", s->id, s->fail->id);
			}
			// Join outputs missing
		}
		list = list_pop ( list );
	}

	list_destroy ( list );
}

// Insert a string to the tree
void ac_addstring ( struct ac_table *g, unsigned int i, unsigned char *string, int m, int alphabet, int *state_transition, unsigned int *state_final ) {

	struct ac_state *state, *next = NULL;
	int j, done = 0;

	// as long as next already exists follow them
	j = 0;
	state = g->zerostate;

	while ( !done && ( next = state->next[*( string + j )] ) != NULL ) {

		state = next;

		if ( j == m )
			done = 1;

		j++;
		
		//printf("character %c state: %i\n", *( string + j ), state->id);
	}

	// not done yet
	if ( !done ) {
		while ( j < m ) {
			// Create new state
			next = malloc ( sizeof ( struct ac_state ) );

			if ( !next )
				fail ( "Could not allocate memory\n" );
				
			next->next = ( struct ac_state ** ) malloc ( alphabet * sizeof ( struct ac_state * ) );

			next->id = g->idcounter++;
			next->output = NULL;
			
			state_transition[state->id * alphabet + *( string + j )] = next->id;
			//printf("setting %i to %i\n", state->id * alphabet + *( string + j ), next->id); 
			
			//printf("Created link from state %i to %i for character %i (j = %i)\n", state->id, next->id, *( string + j ), j );

			//Set all alphabet bytes of the next node's->next to 0
			//This is the _extended_ Aho-Corasick algorithm. A complete automaton is used where all states 
			//have an outgoing transition for every alphabet character of the alphabet
			memset ( next->next, 0, alphabet * sizeof ( struct ac_state * ) );

			state->next[*( string + j )] = next;
			state = next;
			
			//printf("character %c state: %i\n", *( string + j ), state->id);
			j++;
		}
	}
	
	//printf("	Currently at state %i\n", state->id);

	//After finishing with the previous characters of the keyword, add the terminal state if it does not exist
	if ( !state->output ) {
	
		//printf("	For pattern %i added the terminal state %i of %i\n", i, state->id, g->patterncounter);
		state_final[state->id] = 1;

		//allocate memory and copy *string to state->output
		state->output = ( unsigned char * ) malloc ( sizeof ( unsigned char ) * m );
		memcpy ( state->output, string, m );
		
		state->keywordline = g->patterncounter;
		
		g->patterncounter++;
	}
}

unsigned int search_ac ( unsigned char *text, int n, struct ac_table *table ) {

	struct ac_state *head = table->zerostate;
	struct ac_state *r, *s;
	
	int column, matches = 0;
	
	r = head;
	
	for ( column = 0; column < n; column++ ) {
	
		while ( ( s = r->next[*( text + column ) ] ) == NULL )
			r = r->fail;
		r = s;
		
		//printf("column %i r->id = %i\n", column, r->id);
			
		if ( r->output != NULL ) {
			matches++;
			//printf("match at %i for r %i\n", column, r->id);
		}
	}

	return matches;
}

struct ac_table *preproc_ac ( unsigned char **pattern, int m, int p_size, int alphabet, int *state_transition, unsigned int *state_supply, unsigned int *state_final ) {

	unsigned int i;
	
	struct ac_table *table;

	// allocate memory for the table

	table = malloc ( sizeof ( struct ac_table ) );

	if ( !table )
		fail ( "Could not initialize table\n" );

	ac_init ( table, alphabet, state_transition );
	
	for ( i = 0; i < p_size; i++ )
		ac_addstring ( table, i, pattern[i], m, alphabet, state_transition, state_final );
		
	ac_maketree ( table, alphabet, state_supply );
	
	return table;
}

void free_ac ( struct ac_table *table, int alphabet ) {

	ac_destroy ( table, alphabet );

	free ( table );
}
