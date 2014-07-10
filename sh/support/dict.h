// dict.[ch]
// an associative memory
//
// (c) 2005, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_AT_- computer.org
//
// BSD License applies

#ifndef WJDB_DICT_H
#define WJDB_DICT_H

#define DICT_TABLE_LEN 128	///< #elements

/** a dictionary maps a key to a value.
 *
 *  we allocate space for keys, but pass foreign pointers directly
 *  the caller must take care not to free pointers still in the table 
 *
 *  index-based functions start counting at 0 */

/// a (key,value) pair as used in the table
struct dictionary {
	char *key;
	int klen;
	void *value;
};

#define STATIC_DICT(name) static struct dictionary name[DICT_TABLE_LEN];

struct dictionary * dict_create(void);
void dict_clear(struct dictionary *dict, int free_values);
void dict_destroy(struct dictionary *dict, int free_values);

// add / replace / del
int dict_insert(struct dictionary *dict, const char *key, void *value);	
int dict_insert_dup(struct dictionary *dict, const char *key, void *value);	
void dict_replace(struct dictionary *dict, const char *key, void *value);
void dict_del(struct dictionary *dict, const char *key);	
void dict_delex(struct dictionary *dict, const char *key, int n);	
struct dictionary * dict_copy(struct dictionary *dict);

int dict_len(struct dictionary *dict);

// lookup by key / value / index
void * dict_lookup(struct dictionary *dict, const char *key);		
void * dict_lookup_ex(struct dictionary *dict, const char *key, int n);
char * dict_rlookup(struct dictionary *dict, void *data);
void * dict_getnth(struct dictionary *dict, int n);

/// get the next used entry (for internal use only)
static inline int 
__dict_getnext(struct dictionary *dict, int i)
{
	while (++i < DICT_TABLE_LEN)
		if (dict[i].klen)
			return i;
	return -1;
}

/// do something for each filled entry
#define dict_foreach(dict, i, outkey, outval) \
	for ((i) =__dict_getnext(dict, -1); \
	     (i) >= 0 && ((outkey) = dict[i].key) && ((outval) = dict[i].value); \
	     (i) = __dict_getnext(dict, i))

/// retrieve consecutive elements
#define dict_foreach_elem(dict, i, elem) \
	for ((i) = __dict_getnext(dict, -1); \
	     (i) >= 0 && ((elem) = &dict[i]) != NULL; \
	     (i) = __dict_getnext(dict, i))

#endif /* WJDB_DICT_H */

