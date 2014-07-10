/**
 * support/multihash.[ch]
 * a hashtable with lists for each element, to allow unlimited 'siblings'
 *
 * (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
 * email at willem _-at-_ computer DOT org
 *
 * 3-clause BSD applies
 * */

#ifndef SL_SUPP_MHASH
#define SL_SUPP_MHASH

#ifdef __KERNEL__
#include <linux/types.h>
#else
#include <stdint.h>
#endif

#include "hashtable.h"
#include "list.h"

struct multihash {
	// each element is taken as the head of a list
	struct hashtable table;
};

static inline void mhash_add(struct multihash * mh, uint16_t key, void *value)
{
	struct list * list;

	list = hash_lookup(&mh->table, key);
	list = list_append(list, value);
	mh->table.table[key] = list;
}

// remove an entry. 
// or remove all entries with key 'key' by passing NULL as value
static inline void mhash_del(struct multihash * mh, uint16_t key, void *value)
{
	struct list * head, * elem;

	head = hash_lookup(&mh->table, key);
	list_foreach(head, elem)
		if (!value || elem->id == value) {
			head = list_remove(elem);
			mh->table.table[key] = head;
			return;
		}
}

// get all matches for mh(key). the iterator subkey starts at 0
void * mhash_get(struct multihash * mh, uint16_t key, uint16_t subkey);

#endif /* SL_SUPP_MHASH */

