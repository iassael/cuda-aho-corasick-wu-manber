// dict.[ch]
// an associative memory
//
// (c) 2005, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_AT_- computer.org
//
// BSD License applies

#ifdef __KERNEL__
#include <linux/slab.h>
#include <linux/string.h>
#include <linux/module.h>
#else
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#endif

#include "log.h"
#include "macros.h"
#include "dict.h"

#define DICT_BYTELEN (sizeof(struct dictionary) * DICT_TABLE_LEN)

struct dictionary * 
dict_create(void)
{
	struct dictionary *dict;
	
	dict = myalloc(DICT_BYTELEN);
	memset(dict, 0, DICT_BYTELEN);
	return dict;
}

static struct dictionary * 
__dict_find_ex(struct dictionary *dict, const char *key, int n)
{
	int i = -1,  occur = -1;
	
	if (!dict)
		return NULL;

	while (occur < n) {
		if (++i == DICT_TABLE_LEN)
			return NULL;
		if (!dict[i].klen)
			continue;
		if (!memcmp(dict[i].key, key, 
			    max(dict[i].klen, (int) strlen(key)))){
			occur++;
		}
	}
	return &dict[i];
}

void * 
dict_lookup_ex(struct dictionary *dict, const char *key, int n)
{
	struct dictionary * elem = __dict_find_ex(dict, key, n);
	if (elem)
		return elem->value;
	return NULL;
}

void * 
dict_lookup(struct dictionary *dict, const char *key)
{
	return dict_lookup_ex(dict, key, 0);
}

char * 
dict_rlookup(struct dictionary *dict, void *data)
{
	int i = -1;

	while (++i < DICT_TABLE_LEN) {
		if (dict[i].value == data)
			return dict[i].key;
	}
	return NULL;
}

int 
dict_insert_dup(struct dictionary *dict, const char *key, void *value)
{
	int i = 0;

	// skip used items
	while(i < DICT_TABLE_LEN && dict[i].klen)
		i++;
	if (i == DICT_TABLE_LEN) {
		sl_log(LOG_WARN, "exhausted dictionary space");
		return -1;
	}

	// fill item
	dict[i].klen = strlen(key) + 1;
	dict[i].key = myalloc(dict[i].klen);
	memcpy(dict[i].key, key, dict[i].klen);
	dict[i].value = value;
	return i;
}

int 
dict_insert(struct dictionary *dict, const char *key, void *value)
{
	if (!key) {
		sl_log(LOG_WARN, "dict insert NULL key thwarted");
		return -1;
	}

	// there are faster alternatives for duplicate checking
	if (dict_lookup(dict, key)){
		sl_log(LOG_WARN, "dictionary collision on %s", key);
		sl_log(LOG_MSG, key);
		return -1;
	}
	return dict_insert_dup(dict, key, value);
}

void 
dict_replace(struct dictionary *dict, const char *key, void *value) 
{
	struct dictionary *item = __dict_find_ex(dict, key, 0);
	if (item)
		item->value = value;
}

int 
dict_len(struct dictionary *dict)
{
	int i, occur = 0;

	for (i = 0; i < DICT_TABLE_LEN; i++)
		if (dict[i].klen)
			occur++;

	return occur;
}

void * 
dict_getnth(struct dictionary *dict, int n)
{
	int i=-1, j=-1;
	while (++i < DICT_TABLE_LEN)
		if (dict[i].klen && ++j == n)
			return dict[i].value;
	return NULL;
}

void 
dict_delex(struct dictionary *dict, const char *key, int n) 
{
	struct dictionary * entry;
	
	if (!(entry = __dict_find_ex(dict, key, n)))
		return;

	myfree(entry->key);
	entry->klen = 0;
	entry->key = NULL;
	entry->value = NULL;
}

void 
dict_del(struct dictionary *dict, const char *key) 
{
	dict_delex(dict, key, 0);
}

void
dict_clear(struct dictionary *dict, int free_values)
{
	struct dictionary *elem;
	int i;
	
	dict_foreach_elem(dict, i, elem) {
		myfree(elem->key);
		elem->klen = 0;
		if (free_values)
			myfree(elem->value);
	}
}

void 
dict_destroy(struct dictionary *dict, int free_values) 
{
	dict_clear(dict, free_values);
	myfree(dict);
}

struct dictionary * 
dict_copy(struct dictionary *dict)
{
	struct dictionary *new;
	int i;
	char *key;
	void *value;

	if ((new = dict_create()))
		dict_foreach(dict, i, key, value)
			dict_insert(new, key, value);

	return new;
}

#ifdef __KERNEL__
EXPORT_SYMBOL(dict_insert);
EXPORT_SYMBOL(dict_replace);
EXPORT_SYMBOL(dict_lookup);
EXPORT_SYMBOL(dict_rlookup);
EXPORT_SYMBOL(dict_del);
EXPORT_SYMBOL(dict_delex);
EXPORT_SYMBOL(dict_destroy);
#endif


