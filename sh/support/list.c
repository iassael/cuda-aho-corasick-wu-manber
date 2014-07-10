// list.[ch]
// a doubly linked list
//
// (c) 2005, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// LGPL license applies

#include "../support/macros.h"
#include "list.h"

struct list * list_insert_sorted(struct list *start, void *id)
{
	struct list *cur, *new;

	new = list_create(id);
	if (!new)
		return NULL;
	
	if (!start) // start of list: update global startnode
		return new;

	// find our spot in the list
	// exception : test cur (we test cur->next in general)
	if (start->id > id){
		start->prev = new;
		new->next = start;
		return new;
	}
	if (start->id == id){
//		printf("skipping duplicate : cur=%p\n",new->id);
		free(new);
		return start;
	}

	cur = start;
	while (cur->next && cur->next->id < id)
		cur = cur->next;
	
	if (!cur->next){ // end of list: append or place just before the end-node
//		printf("inserting (%p) : cur=%p new=%p\n",id,cur ? cur->id : "[ ]",new->id);
		cur->next = new;
		new->prev = cur;
		return start;
	}
	
	if (cur->next->id == id){ // exception : found a duplicate. remove
//		printf("skipping duplicate : cur=%p\n",new->id);
		free(new);
		return start;
	}
	
	
	// insert into the sorted list
//	printf("inserting (%p) : cur=%p new=%p cur->next=%p\n",id,cur ? cur->id : "[ ]",new->id, cur->next->id);
	new->next = cur->next;
	new->prev = cur;
	cur->next->prev = new;
	cur->next = new;
	
	return start;
}

