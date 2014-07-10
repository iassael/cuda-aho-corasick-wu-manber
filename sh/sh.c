#include "../smatcher.h"

/// free an AC table from a given startnode (recursively)
void sh_free ( struct ac_state *state, int alphabet ) {
	
	int i;

	for ( i = 0; i < alphabet; i++ )
		if ( state->next[i] )
			sh_free ( state->next[i], alphabet );

	if ( state->output )
		free ( state->output );
	
	free ( state->next );
	free ( state );
}

/// initialize the empty-table
void sh_init ( struct ac_table *g, int alphabet, int *state_transition ) {

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
void sh_destroy ( struct ac_table *in, int alphabet ) {

	int i;

	for ( i = 0; i < alphabet; i++ )
		if ( in->zerostate->next[i] && in->zerostate->next[i]->id > 0 ) {
			
			//printf("id: %i i: %i\n", in->zerostate->next[i]->id, i);
		
			sh_free ( in->zerostate->next[i], alphabet );
			in->zerostate->next[i] = NULL;
		}
	free ( in->zerostate->next );
	free ( in->zerostate );
}

// Insert a string to the tree
void sh_addstring ( struct ac_table *g, unsigned int i, unsigned char *string, int m, int alphabet, int *state_transition, unsigned int *state_final ) {

	struct ac_state *state, *next = NULL;
	int j, done = 0;

	// as long as next already exists follow them
	j = m - 1;
	state = g->zerostate;

	while ( !done && ( next = state->next[*( string + j )] ) != NULL ) {
	
		//printf("id: %i j: %i\n", state->id, j);

		state = next;

		if ( j <= 0 )
			done = 1;

		j--;
		
		//printf("character %c state: %i\n", *( string + j ), state->id);
	}

	// not done yet
	if ( !done ) {
		while ( j >= 0 ) {
			// Create new state
			next = malloc ( sizeof ( struct ac_state ) );

			if ( !next )
				fail ( "Could not allocate memory\n" );
				
			next->next = ( struct ac_state ** ) malloc ( alphabet * sizeof ( struct ac_state * ) );

			next->id = g->idcounter++;
			next->output = NULL;
			
			state_transition[state->id * alphabet + *( string + j )] = next->id;
			
			//printf("Created link from state %i to %i for character %c  (j = %i)\n", state->id, next->id, *( string + j ), j );

			//Set all alphabet bytes of the next node's->next to 0
			//This is the _extended_ Aho-Corasick algorithm. A complete automaton is used where all states 
			//have an outgoing transition for every alphabet character of the alphabet
			memset ( next->next, 0, alphabet * sizeof ( struct ac_state * ) );

			state->next[*( string + j )] = next;
			state = next;
			
			j--;
		}
	}

	//After finishing with the previous characters of the keyword, add the terminal state if it does not exist
	if ( !state->output ) {
	
		state_final[state->id] = 1;

		//allocate memory and copy *string to state->output
		state->output = ( unsigned char * ) malloc ( sizeof ( unsigned char ) * m );
		memcpy ( state->output, string, m );
		
		//printf("Adding output %s to state %i\n", state->output, state->id);

		state->keywordline = g->patterncounter;
		
		g->patterncounter++;
	}
}

unsigned int search_sh ( int m, unsigned char *text, int n, struct ac_table *table, int *bmBc ) {

	struct ac_state *head = table->zerostate;
	struct ac_state *r, *s;
	
	int column = m - 1, matches = 0, j;
	
	r = head;

	while ( column < n ) {
		
		r = head;
		j = 0;

		while ( j < m && ( s = r->next[*( text + column - j )] ) != NULL ) {
				
			r = s;
			j++;
		}
			
		if ( r->output != NULL )
			matches++;
			
		column += bmBc[text[column]];
	}
	
	return matches;
}

struct ac_table *preproc_sh ( unsigned char **pattern, int m, int p_size, int alphabet, int *state_transition, unsigned int *state_final ) {

	unsigned int i;
	
	struct ac_table *table;

	// allocate memory for the table

	table = malloc ( sizeof ( struct ac_table ) );

	if ( !table )
		fail ( "Could not initialize table\n" );

	sh_init ( table, alphabet, state_transition );

	for ( i = 0; i < p_size; i++ )
		sh_addstring ( table, i, pattern[i], m, alphabet, state_transition, state_final );

	return table;
}

void free_sh ( struct ac_table *table, int alphabet ) {

	sh_destroy ( table, alphabet );

	free ( table );
}

