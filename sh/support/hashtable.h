// hashtable.[ch]
// a double-hashed hashtable
//
// (c) 2005, willem de bruijn, vrije universiteit amsterdam
// email at willem -_at_- computer.org
//
// BSD license applies

#ifndef WJDB_HASHTABLE
#define WJDB_HASHTABLE

#include "macros.h"

#define HASHTBL_LEN 	431	// use a prime
#define MAX_DOUBLEHASH	24	// stop searching after so many rounds
struct hashtable {
	// TODO: add length field, make default size (much) smaller and enable resizing
	void * table[HASHTBL_LEN];
};

#define hash_lookup_fast(htable, key) \
	((htable)->table[key])

// check against overflow. constructed in such a way that
// we can see in what function the overflow occurred.
#ifndef NDEBUG
#define hash_lookup(htable, key) \
	(( ((unsigned long) key) < HASHTBL_LEN) ? \
	       hash_lookup_fast(htable, key) : \
	       ((void*) (dprintf("BUG: key %d out of bounds in %s\n", \
			 key, __FUNCTION__) & 0L)))

#define hash_insert_at_unconditional(htable, value, key) \
	do {(htable)->table[key] = value;} while (0)

// insert and check against overwriting. also see hash_lookup
// returns the key, or <0 on error
#define hash_insert_at(htable, value, key) \
	((!(htable)->table[key]) ? \
	 	(((htable)->table[key] = value) ? key : -1) : \
		((dprintf("BUG: key %d in use in %s.%d\n", \
			 key, __FUNCTION__, __LINE__) & 0L)))
#else

#define hash_lookup(htable, key) \
	((((unsigned long) key) < HASHTBL_LEN) ? \
	       hash_lookup_fast(htable, key) : 0)

// insert and check against overwriting. also see hash_lookup
#define hash_insert_at(htable, value, key) \
	((!(htable)->table[key]) ? \
	 	(((htable)->table[key] = value) ? key : -1) : 0)
#endif

static inline int
hash_calc(void * value, int runno, int maxhash)
{
	int h, k, i;
	
	i = 0;
	h = ((unsigned long) value) % maxhash;		// primary hash function
	k = ((unsigned long) value) % (maxhash - 2);	// secondary hash function
	
	return (h + runno * k) % maxhash;
}

int hash_lookup_by_value(struct hashtable *hash, void * value, int nth);

// the __FUNCTION__ helps me locate collision origins
#define hash_insert(a,b) __hash_insert(a,b,__FUNCTION__)
int __hash_insert(struct hashtable *hash, void * value, const char *func);

static inline int 
hash_del(struct hashtable *hash, int key)
{
#ifndef NDEBUG
	check(key >= 0 && key < HASHTBL_LEN);
#endif
	hash->table[key] = NULL;
	return 0;
}

/** use the hashtable as a simple list */
static inline 
int hash_getnext(struct hashtable *hash, int key)
{
	if (key >= HASHTBL_LEN || key < -1)
		return -1;
	
	while (!hash->table[++key])
		if (key == HASHTBL_LEN-1)
			return -1;

	return key;
}

// use an integer for key
#define hash_foreach(table, key, ptr) \
	for(key = hash_getnext(table, -1);\
	    key >= 0 && key < HASHTBL_LEN && \
	    (((ptr) = hash_lookup_fast(table,key)) != NULL);\
	    key = hash_getnext(table,key))

// return 0 if the hashtable contains a value, !0 (i.e., true) otherwise
static inline int 
hash_empty(struct hashtable *table)
{
	int i;

	i = hash_getnext(table, -1);
	return (i < 0);
}

static inline int 
hash_len(struct hashtable *table)
{
	int i, count=0;
	
	i = hash_getnext(table,-1);
	while(i >= 0){
		count++;
		i = hash_getnext(table,i);
	}

	return count;
}

#endif /* WJDB_HASHTABLE */

