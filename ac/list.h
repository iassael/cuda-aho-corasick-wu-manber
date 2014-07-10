// list.[ch]
// a doubly linked list
//
// (c) 2005, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// BSD license applies

#ifndef WDB_SLIST_H
#define WDB_SLIST_H

#ifdef __KERNEL__
#include <linux/kernel.h>
#include <linux/slab.h>
#include <linux/string.h>
#else
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#endif


struct list {
	void *id;
	struct list *next;
	struct list *prev;
};

static inline struct list * list_create(void *id)
{
	struct list * new;

	new = malloc (sizeof(struct list));
	if (!new)
		return NULL;
	new->id = id;
	new->next = NULL;
	new->prev = NULL;

	return new;
}

static inline struct list * list_insert(struct list *start, void *id)
{
	struct list *new = list_create(id);

	if (!new)
		return NULL;
	
	if (start){
		new->next = start;
		start->prev = new;
	}
	return new;
}

static inline struct list * list_append(struct list *start, void *id)
{
	struct list *new = list_create(id);
	struct list* cur;

	if (!new)
		return NULL;
	
	if (!start)
		return new;
	
	cur = start;
	while (cur->next)
		cur = cur->next;
	cur->next = new;
	new->prev = cur;
	return start;
}

/** strange function for a list 
  * used only for duplicate removal 
  * 
  * note that the function returns NULL in two
  * distinct cases: no 'start', or 'start' is the only item
  */
static inline struct list * list_pop(struct list *start)
{
	struct list *tmp;

	if (!start)
		return NULL;
	
	tmp = start;
	start = start->next;
	free(tmp);
	
	return start;
}

static inline struct list * list_invert(struct list *start)
{
	struct list *cur, *tmp=NULL;
	
	if (!start->next)
		return start;

	cur = start;
	// swap {prev,next} pointers
	while (cur){
		tmp = cur->next;
		cur->next = cur->prev;
		cur->prev = tmp;
		tmp = cur;
		cur = cur->prev;
	}

	return tmp;
}

// return the item in the list that matches the id
static inline struct list * list_exists(struct list *start, void * id)
{
	struct list *cur;

	if (!start)
		return NULL;
	
	// find our spot in the list
	cur = start;
	while (cur && cur->id != id)
		cur = cur->next;
	if (!cur)
		return NULL;
	else	
		return cur;
}

/** unlink an item. can be used together with list_foreach */
static inline struct list * list_unlink(struct list *cur)
{
	struct list *tmp = NULL;

	if (cur->next){
		cur->next->prev = cur->prev;
		tmp = cur->next;
	}
	if (cur->prev){
		cur->prev->next = cur->next;
		tmp = cur->prev;
	}

	if (!tmp)
		return NULL; // no cur->next && no cur->prev ? then it's an empty list
		
	while (tmp->prev)
		tmp = tmp->prev;
	return tmp; // return the new startnode
}

/** remove id if it exists. returns start of the list */
static inline struct list * list_remove(struct list *cur)
{
	struct list * elem = list_unlink(cur);
	free(cur);
	return elem;
}

static inline struct list * list_remove_id(struct list * list, void * id) 
{
	struct list * elem = list_exists(list, id);
	if (elem)
		return list_remove(elem);
	return list;
}

struct list * list_insert_sorted(struct list *start, void *id);

static inline int list_len(struct list *list)
{
	int i=0;
	while (list){
		i++;
		list = list->next;
	}
	return i;
}

#define list_foreach(list, cur) \
	for (cur = list; cur; cur = cur->next)

#define list_destroy(deadlist) \
	while (deadlist) deadlist = list_pop(deadlist)

#endif /* WDB_SLIST_H */

