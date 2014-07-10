#include "../smatcher.h"
/*
//Create the first node of the list
struct node* setup_head ( char label ) {
	
	struct node* newState = malloc ( sizeof ( struct node ) );
	
	newState->label = label;
	newState->id = 0;
	newState->supply = NULL;
	newState->next = NULL;

	return newState;
}

//append a node and visit it
void append_node ( struct node** lastState, char label, int id ) {
	
	struct node* newState = malloc ( sizeof ( struct node ) );

	newState->label = label;
	newState->id = id;
	newState->supply = NULL;
	newState->next = NULL;
	
	//the ->next of the last state in the list will now point to this state
	(*lastState)->next = newState;
	
	//The current node in the list will now be the last node we created
	*lastState = newState;
}

void free_kmp ( struct node* state ) {

	struct node* tmp;
	
	while ( 1 ) {
	
		tmp = state;
		state = state->next;
		
		free ( tmp );
	
		if ( state == NULL )
			break;
	}
}

//create the supply link for a node
void addSupply ( struct node* head, int current, struct node* supply ) {

	int i;

	struct node* currentState = head;
	
	for ( i = 0; i < current; i++ )
		currentState = currentState->next;
		
	currentState->supply = supply;

}

struct node* preKmpList ( struct node* head, unsigned int *pattern, int m ) {

	int i = 0, k;
	
	struct node* currentState = head;
	
	struct node* j = NULL;
	
	for ( k = 1; k <= m; k++ )
		append_node( &currentState, pattern[k], k );
		
	while (i < m) {
	
		while ( j != NULL && pattern[i] != j->label)
			j = j->supply;
		
		i++;

		if ( j == NULL )
			j = head;
		else
			j = j->next;
			
		if ( i < m && pattern[i] == j->label )
			addSupply( head, i, j->supply );
		else	
			addSupply( head, i, j );
	}
	
	return head;
}

unsigned int searchList ( struct node* head, unsigned int *pattern, int m, unsigned int *text, int n ) {

	int i = 0;

	struct node* j = head;

	while (i < n) {

		//mismatch occurs
		while ( j != NULL && j->label != text[i] )
			j = j->supply;

		i++;
		
		if ( j == NULL )
			j = head;
		else
			j = j->next;
		
		if ( j->id >= m ) {
			return ( i - j->id );
			
			printf("->%i\n", i - j->id);

			j = j->supply;
		}
	}
}
*/

void preKmp ( int *next, unsigned char *p, int m ) {

	int i=0;
	int j=-1;
	next[0] = -1;

	while (i < m) {
		
		while ( j >= 0 && p[i]!=p[j] )
			j = next[j];

		i++; j++;
		
		if ( i < m && p[i] == p[j] )
			next[i] = next[j];
		else
			next[i] = j;
	}
}
/*
void search ( int *next, unsigned char *pattern, int m, unsigned char *text, int n ) {

	int i = 0;
	int j = 0;

	while (i < n) {

		//mismatch occurs
		while (j >= 0 && pattern[j] != text[i])			
			j = next[j];
			
		i++;
		j++;

		if (j >= m) {
			//printf("->%i\n", i - j);
			j = next[j];
		}
	}
}

int main ( void ) {

	int i;

	int m = 8;
	unsigned char *pattern = (unsigned char *)"AACGTAAC";
	
	int n = 12;
	unsigned char *text = (unsigned char *)"TAATAACGTAAC";
	
	preKmp( pattern, m );
	
	search( pattern, m, text, n );
	
	for ( i = 0; i < m; i++ )
		printf("%i\n", next[i]);
		
	printf("\n");
	
	struct node* head = setup_head( pattern[0] );
	
	struct node* state = head;
		
	preKmpList( state, pattern, m );
	
	searchList ( state, pattern, m, text, n );
	
	while ( 1 ) {

		if ( state->label ) {
			if ( state->supply != NULL )
				printf("Node %c points to node %i with a label %c\n", state->label, state->supply->id, state->supply->label);
			else
				printf("Node %c points to NULL\n", state->label);
		}	
		state = state->next;
		
		if ( state == NULL )
			break;
	}
	
	free_kmp ( head );
	
	return 0;
}
*/

