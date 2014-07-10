/* support/slist.[ch]
 * a stack-based list implementation
 * 
 * (c) 2008, Willem de Bruijn, Vrije Universiteit Amsterdam
 * GPLv2 applies
 * 
 * */

#ifdef __KERNEL__
#include <linux/kernel.h>
#include <linux/string.h>
#else
#include <stdlib.h>
#include <string.h>
#endif

#include "slist.h"

#define INCFACTOR 2	/**< expansion rate */
#define STARTLEN  4

int 
slist_add(struct slist * sl, unsigned long key, void *arg)
{
	// realloc
	if (sl->used == sl->len) {
		struct slist_elem *bak;
		int bytelen;

		bytelen = sl->len * sizeof(struct slist_elem);
		if (bytelen) {
			bak = sl->list;
			sl->list = myalloc(INCFACTOR * bytelen);
			memcpy(sl->list, bak, bytelen);
			myfree(bak);
			sl->len *= INCFACTOR;
		}
		else {
			sl->len = STARTLEN;
			sl->list = myalloc(sl->len * sizeof(struct slist_elem));
		}
	}

	// add
	sl->list[sl->used].key = key;
	sl->list[sl->used].arg = arg;
	sl->used++;
	return 0;
}

int
slist_del(struct slist *sl, unsigned long key)
{
	int i = __sllist_get(sl, key);

	if (i < 0) {
		sl_log(LOG_WARN, "deallocation from slist failed");
		return -1;
	}

	// place last element into newly created hole
	if (i < sl->used - 1) {
		sl->list[i].key = sl->list[sl->used - 1].key;
		sl->list[i].arg = sl->list[sl->used - 1].arg;
	}
	sl->used--;
	
	return 0;
}

