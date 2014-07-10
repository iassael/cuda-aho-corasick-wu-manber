// hashtable.[ch]
// a double-hashed hashtable
//
// (c) 2008, willem de bruijn, vrije universiteit amsterdam
// email at willem -_at_- computer.org
//
// BSD license applies

#include "hashtable.h"

/** lookup an element by calling the hashfunction, like in hash_insert. 
  return the nth match. NB: n starts at 1! */
int 
hash_lookup_by_value(struct hashtable *hash, void * value, int nth) 
{
	int i=-1, match=0, key=-1;
	
	check(nth > 0);
	while (i < MAX_DOUBLEHASH && match < nth){
		key = hash_calc(value, ++i, HASHTBL_LEN);
		if (hash->table[key] == value)
			match++;
	}
	if (match == nth)
		return key;
	else
		return -1;
}

/** insert an item, we use double hashing for collision resolution */
int 
__hash_insert(struct hashtable *hash, void * value, const char *func)
{
	int i=0, index;
	
	if (!value)
		returnbug(-1);
		
	index = hash_calc(value, i, HASHTBL_LEN);
	while (hash->table[index] && i < MAX_DOUBLEHASH) {
		index = hash_calc(value, ++i, HASHTBL_LEN);
#ifndef NDEBUG
		if (hash->table[index] == value)
			dprintf("warning : duplicate hash %d: %p==%ld in %s\n",
				index, value, (long) value, func);		
#endif	
	}
		
	if (i == MAX_DOUBLEHASH) { // give up
		dprintf("(BUG) hash full in %s\n", func);
		return -1;
	}
	
	hash->table[index] = value;
	return index;
}

