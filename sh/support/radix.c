/* support/radix.[ch] 
 * implementation of a radix tree
 * the tree works on arbitrary binary strings. \0 is not necessary
 * 
 * (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
 * GPLv2 applies
 *
 * TODO: replace per-node allocation with page-sized memory pools
 
 */

#ifdef __KERNEL__
#include <linux/slab.h>
#include <linux/string.h>
#else
#include <stdlib.h>
#include <string.h>
#endif

#include "macros.h"
#include "radix.h"

// TODO: allow bitlengths != 8
#define BITLEN 8		   // number of bits per character
#define ALPHABET_LEN (1 << BITLEN) // number of unique letters in the alphabet

struct radix_node {
	struct radix_node *parent; 
	struct radix_node *children[ALPHABET_LEN];

	void *value;
	
	// radix-tree specific: allow path compression in "remainder"
	unsigned char *remainder;
	size_t rlen;
};

// return the node furthest on the path to key
// if on return klen > 0 this is the remainder that we failed to resolve
// if it is 0 the returned node is the node associated with key
static struct radix_node *
__radix_lookup(struct radix_node *root, unsigned char *key, size_t *klen)
{
	int koff = 0, len = *klen;

	if (!root)
		return NULL;

	// traverse the tree down as far as possible
	while (koff < len && root->children[(int) key[koff]]) {
		root = root->children[(int) key[koff]];
		koff++;
	}
	len -= koff;

	// is node an internal node or a compressed leaf?
	if (len) {
		// if so, is it the compressed leaf we want?
		if (root->rlen == len &&
		    !memcmp(root->remainder, key + koff, len))
			len = 0;
	}
	
	*klen = len; // ILLEGAL! TODO: FIX
	return root;
}

/// return the value associated with key, if any
void * radix_lookup(struct radix_node *root, unsigned char *key, size_t klen)
{
	root = __radix_lookup(root, key, &klen);
	if (!root || klen)
		return NULL;
	return root->value;
}

/// return the rightmost child of node, if any
static struct radix_node * __radix_max_child(struct radix_node *node)
{
	int letter = ALPHABET_LEN;
	while (--letter >= 0 && node->children[letter] == NULL) {}
	
	if (letter < 0)
		return NULL;
	return node->children[letter];
}

// find the rightmost leaf of node, if any
static struct radix_node * __radix_max_descendant(struct radix_node *node)
{
	struct radix_node *child;

	// if succesfull, recurse down to the rightmost child
	while ((child = __radix_max_child(node)) != NULL)
		node = child;
	return node;
}

/// find the largest leaf smaller than node (which may be an internal node)
static struct radix_node *
__radix_lookup_predecessor(struct radix_node *node)
{
	int letter;

	do {
		if (!node->parent)
			return NULL;

		// first find 'node' again in the list
		letter = ALPHABET_LEN - 1;
		while (node->parent->children[letter] != node)
			letter--;
		
		// find the largest smallest letter of the same parent
		while (letter > 0) {
			letter--;
			if (node->parent->children[letter])
				return __radix_max_descendant(
						node->parent->children[letter]);
		};

		// not found? move up the tree
		// special case the root node
		if (!node->parent)
			return __radix_max_descendant(node);
		node = node->parent;
	} while (node);

	return NULL; // never reached
}

/// return the value associated with key or, if there is no perfect match,
/// the nearest predecessor
void * radix_lookup_predecessor(struct radix_node *root, unsigned char *key, 
				size_t klen)
{
	struct radix_node * pre;
	size_t rest = klen;
	int letter;

	root = __radix_lookup(root, key, &rest);
	if (!root || !rest) // error or perfect match
		return root;

	// find the largest child smaller than us (at the current level)
	letter = key[klen - rest];
	while (letter >= 0) {
		pre = root->children[letter--];
		if (pre) { // found
			return __radix_max_descendant(root)->value;
		}
	}

	// no smaller element of the same parent? try parent's predecessor
	pre = __radix_lookup_predecessor(root);
	if (pre)
		return pre->value;
	else
		return NULL;
}

static struct radix_node * __radix_alloc(void)
{
	struct radix_node *child;

	child = myalloc(sizeof(struct radix_node));
	memset(child, 0, sizeof(struct radix_node));
	return child;
}

/// add an intermediate node
static struct radix_node * 
__radix_add_intermediate(struct radix_node *parent, unsigned char key)
{
	struct radix_node *child = __radix_alloc();

	// connect it to the parent
	parent->children[key] = child;
	child->parent = parent;
	return child;
}

/// add a rootnode : a leaf without a parent
static struct radix_node * 
__radix_add_root(unsigned char *key, size_t klen, void *value)
{
	struct radix_node * root = __radix_alloc();

	root->rlen = klen;

	root->remainder = myalloc(root->rlen + 1);
	memcpy(root->remainder, key, root->rlen);
	root->remainder[root->rlen] = '\0';
	
	root->value = value;
	return root;
}

/// add a node with a completely compressed path to a node
static struct radix_node * 
__radix_add_leaf(struct radix_node *parent, unsigned char *key, size_t klen, 
		 void *value)
{
	struct radix_node *child;

	// create the child
	// if the length > 1 use a compressed node, else don't
	if (klen > 1) {
		child = __radix_add_root(key + 1, klen - 1, value);
		parent->children[key[0]] = child;
		child->parent = parent;
	}
	else {
		child = __radix_add_intermediate(parent, key[0]);
		child->value = value;
	}

	return child;
}

/// return the new node on success, NULL on failure
struct radix_node * radix_insert(struct radix_node *root, unsigned char *key, 
				 size_t klen, void *value)
{
	struct radix_node *node;
	int koff = 0;

	if (!klen)
		return NULL;

	// special case: no root
	if (!root)
		return __radix_add_root(key, klen, value);	
		
	// descend down the tree as far as possible towards our goal
	while (koff < klen && root->children[key[koff]]) {
		root = root->children[key[koff]];
		koff++;
	}

	// add the new compressed node and any necessary intermediates 
	
	// split existing compressed node in shared subpath and 
	// private smaller compressed node
	if (root->rlen) {
		int i = 0, k;

		// see to what extend the two substrings overlap
		while (i < root->rlen 
		       && koff < klen
		       && root->remainder[i] == key[koff]) {
			i++;
			koff++;
		}
		
		// special case 1: new node is identical to existing: fail
		if (koff == root->rlen)
			return NULL;

		// add intermediate nodes
		k = 0;
		node = root;
		while (k < i)
			node = __radix_add_intermediate(node, 
							root->remainder[k++]);
		
		// add new leaf for pre-existing value
		__radix_add_leaf(node, root->remainder + i, 
			         root->rlen - i, root->value);
	
		// free stale compressed data
		myfree(root->remainder);
		root->rlen = 0;
		root = node;
	}
	
	// add new value
	// should we add a new leaf, or is this an internal node?
	if (koff == klen) {
		// disallow duplicates
		if (root->value)
			return NULL;
		else {
			root->value = value;
			return root;
		}
	}
	else
		return __radix_add_leaf(root, key + koff, klen - koff, value);
}

// TODO : recompress nodes that no longer have any siblings
int radix_delete(struct radix_node *root, unsigned char *key, size_t klen)
{
	struct radix_node *node, *parent;
	int letter = ALPHABET_LEN - 1;

	node = __radix_lookup(root, key, &klen);
	if (!node || klen)
		return 0;
	
	// free the node and its subtree
	parent = node->parent;
	radix_destroy(node);

	// have we removed the last element?
	if (!parent)
		return 1;
		
	// remove the link in the parent node 
	while (parent->children[letter] != node)
		letter--;
	parent->children[letter] = NULL;
	return 0;
}

void radix_destroy(struct radix_node *root)
{
	int letter;

	// first recurse into subtrees
	for (letter = 0; letter < ALPHABET_LEN; letter++)
		if (root->children[letter])
			radix_destroy(root->children[letter]);

	// then destroy the node itself
	if (root->rlen)
		myfree(root->remainder);
	myfree(root);
}


