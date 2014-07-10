/**
 * support/multihash.[ch]
 * a hashtable with lists for each element, to allow unlimited 'siblings'
 *
 * (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
 * email at willem _-at-_ computer DOT org
 *
 * 3-clause BSD applies
 * */

#include "multihash.h"

// mhash_get without optimizations
static inline struct list * __mhash_getelem(struct multihash * mh, 
			       	            uint16_t key, 
	      	 		            uint16_t subkey)
{
	struct list * head, * elem;
	
	head = hash_lookup(&mh->table, key);
	list_foreach(head, elem)
		if (!subkey--)
			return elem;
	return NULL;
}

// optimization: expect this request to be part of a loop
// we cache the pointer and compare {mh, key, subkey} with cached version
//
// nb: this is unsafe in a very specific situation: when the list into
// which the cached pointer points is changed between calls. This is
// highly unlikely. TODO: make certain this cannot occur
void * mhash_get(struct multihash * mh, uint16_t key, uint16_t subkey)
{
	static struct multihash * s_mh;
	static uint16_t s_key;
	static uint16_t s_sub;
	static struct list *s_elem;
	struct list *elem;

	// next iterator element in the current loop?
	if (s_mh == mh && s_key == key && s_sub + 1  == subkey) {
		s_sub++;
		s_elem = s_elem->next;
		elem = s_elem;
	} 
	else if (!subkey) { // start of a new loop?
		s_mh = mh;
		s_key = key;
		s_sub = 0;
		s_elem = __mhash_getelem(mh, key, 0);
		elem = s_elem;
	} else {
		elem = __mhash_getelem(mh, key, subkey);
	}
	
	if (elem)
		return elem->id;
	else
		return NULL;
}

#ifdef __KERNEL__

#include <linux/module.h>
EXPORT_SYMBOL(mhash_get);

#endif

