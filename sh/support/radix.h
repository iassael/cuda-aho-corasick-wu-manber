/* support/radix.[ch] 
 * implementation of a radix tree
 *
 * the tree works on arbitrary binary strings. \0 is not necessary
 * duplicate keys are not allowed
 * 
 * (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
 * GPLv2 applies
 */

#ifndef SL_SUPPORT_RADIX_H
#define SL_SUPPORT_RADIX_H

#ifdef __KERNEL__
#include <linux/types.h>
#else
#include <sys/types.h>
#endif

struct radix_node;

/// lookup a value by a <key, klen> tuple. 
//  returns NULL on failure
void * radix_lookup(struct radix_node *, unsigned char *, size_t);

/// lookup a value or its nearest predecessor.
//  returns NULL on failure
void * radix_lookup_predecessor(struct radix_node *, unsigned char *, size_t);

/// insert a <key, klen, value> tuple. 
//  returns NULL on failure
struct radix_node * radix_insert(struct radix_node *, unsigned char *, size_t, 
				 void *);

/// delete a <key, klen, value> tuple by passing the associated node
//  returns 1 if we removed the rootnode, 0 otherwise
int radix_delete(struct radix_node *, unsigned char *key, size_t keylen);

/// destroy an entire tree
void radix_destroy(struct radix_node *);

#endif

