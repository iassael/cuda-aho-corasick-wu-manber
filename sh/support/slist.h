/* support/slist.[ch]
 * A stack-based list implementation optimized for lookup
 * (at the cost of insertion and deletion).
 * 
 * (c) 2008, Willem de Bruijn, Vrije Universiteit Amsterdam
 * GPLv2 applies
 * 
 * */

#include "macros.h"
#include "log.h"

struct slist_elem {
	unsigned long key;
	void *arg;
};


/* A stack-based list allocates an array of pointers
 * and grows as needed. 
 *
 * No initialization is necessary besides setting len and used to 0.
 * */
struct slist {
	int len;
	int used;

	struct slist_elem *list;
};

static inline int
__sllist_get(struct slist *sl, unsigned long key)
{
	int i;

	for (i = 0; i < sl->used; i++) {
		if (key == sl->list[i].key)
			return i;
	}

	return -1;
}

static inline void *
slist_get(struct slist * sl, unsigned long key)
{
	int i = __sllist_get(sl, key);

	if (likely(i >= 0))
		return sl->list[i].arg;
	else
		return NULL;
}

int slist_add(struct slist * sl, unsigned long key, void *elem);
int slist_del(struct slist *sl, unsigned long key);

