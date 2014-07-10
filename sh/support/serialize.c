// serialize.[ch]
// pack/unpack a bunch of parameters
//
// (c) 2005, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// BSD license applies

#ifndef __KERNEL__
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <sys/types.h>
#else
#include <linux/kernel.h>
#include <linux/types.h>
#include <linux/string.h>
#include <linux/slab.h>
#include <linux/module.h>
#endif

#include "macros.h"
#include "serialize.h"

//a maximum, this (theoretically) limits stringlen, and can serve as a rudimentary test
#define MAXTUPLELEN 1 << 16

inline int serialized_totallen(const char *serialized)
{
	struct shead *shead;
    	
	if (!serialized)
	    return 0;
	shead = (struct shead*) serialized;
	return sizeof(struct shead) + shead->bytelen + (sizeof(uint32_t) * shead->tuplelen);
}

char * serialize(int argcount, ...)
{
	__attribute__((unused)) va_list list;
	int i, empty;
	uint32_t elemsize, offset; 
	char *output, *param;
	struct shead shead;

	// phase 1 : figure out how much space we need to allocate
	shead.bytelen = 0;
	shead.tuplelen = 0;
	va_start (list, argcount);
	for (i = 0; i < (2 * argcount); i += 2) {
		shead.bytelen += va_arg(list, uint32_t);
		shead.tuplelen++;
		empty = va_arg(list, int); // skip the data. gives warning with GCC 4.1. FIXED
	}
	va_end (list);	
		
	// phase 2 : copy the data
	output = myalloc (sizeof(struct shead) + (argcount * sizeof(uint32_t)) + shead.bytelen);
	if (!output){
		dprintf("failed to allocated memory (size is %d)\n",shead.bytelen);
		return NULL;
	}

	// write the string header
	offset = sizeof(struct shead);
	memcpy(output, &shead, offset);
	
	// place the elements
	va_start (list, argcount);
	for (i=0;i<argcount;i++){
		// first write the header. this is always 32 bits long
		elemsize = (uint32_t) va_arg(list, int);
		memcpy(output + offset, &elemsize, sizeof(uint32_t));
		offset += sizeof(uint32_t);
		// now write the data
		param = va_arg(list, char*);
		if (elemsize){
			memcpy(output + offset, param, elemsize);
			offset += elemsize; // offset into the output buffer
		}
	}
	va_end (list);

	return output;
}

char * serialize_add(char *serialized, int dlen, const void *data)
{
	struct shead *shead;
	char *new;
	int oldlen;
	uint32_t newheader;
	
	shead = (struct shead *) serialized;
	if (!serialized || shead->tuplelen > MAXTUPLELEN) // some safety integrity checks. 
		return NULL;
	
	oldlen = serialized_totallen(serialized);
	newheader = (uint32_t) dlen;
	new = myalloc (oldlen + dlen + sizeof(uint32_t));
	if (!new)
		return NULL;

	memcpy(new,serialized,oldlen);
	memcpy(&new[oldlen],&newheader,sizeof(uint32_t));
	memcpy(&new[oldlen+sizeof(uint32_t)],data,dlen);
	shead = (struct shead *) new;
	shead->tuplelen++;
	shead->bytelen+=dlen;
	myfree(serialized);
	return new;
}

int is_serialized(const char *serialized)
{
	struct shead *shead = (struct shead *) serialized;
	if (!shead)
		return 0;

	if (shead->tuplelen < MAXTUPLELEN && shead->tuplelen <= shead->bytelen)
		return 1;
	else
		return 0;
}

int deserialize(char *string, ...)
{
	va_list list;
	struct shead shead;
	int i, j, offset;
	uint32_t elemsize;
	char **ptr;

	memcpy(&shead, string, sizeof(struct shead));
	offset = sizeof(struct shead);
	va_start (list, string);
	for (i=0;i<shead.tuplelen;i++){
		memcpy(&elemsize, string + offset, sizeof(uint32_t));
		offset += sizeof(uint32_t);
		ptr = va_arg(list, char**);
		*ptr = myalloc(elemsize);
		if (!*ptr)
			goto cleanup;
		memcpy(*ptr, string + offset, elemsize);
		offset += elemsize;
	}
	va_end (list);

	return 0;

cleanup:
	// free previously allocated members
	va_end (list);
	va_start (list, string);
	for(j=0;j<i;j++){
		ptr = va_arg(list, char**);
		myfree (*ptr);
	}
	return -1;
}

const char * serialized_data(const char *serialized, int elemno, int *dlen)
{
	int i;
	uint32_t offset;
	uint32_t itemlen;
	struct shead *shead;
	
	shead = (struct shead *) serialized;
	check_ptr(shead->tuplelen > elemno);
		
	offset = sizeof(struct shead);
	for(i = 0; i < elemno; i++){
		memcpy(&itemlen, serialized + offset, sizeof(uint32_t));
		offset += sizeof(uint32_t) + itemlen;
	}
	if (dlen)
		// BUG: uint32_t -> int
		memcpy(dlen, serialized + offset, sizeof(uint32_t));
	if (*(uint32_t*) serialized + offset == 0) // NULL pointer?
		return NULL;
	return &serialized[offset + sizeof(uint32_t)]; // skip past header
}

char * serialize_duplicate(char * in)
{
    char *out;
    int size;
    
    check_ptr (in && is_serialized(in));

    size = serialized_totallen(in);
    out = myalloc (size);
    check_ptr(out);

    memcpy(out, in, size);
    return out;
}

char * serialize_merge(char * one, char * two, int del)
{
    struct shead *h_out, *h_one, *h_two;
    char * out;
    int size_one, size_two;
    
    if (!one && !two)
		return NULL;
    if (!one){
		if (del)
			return two;
		else
			return serialize_duplicate(two);
    }
    if (!two){
		if (del)
			return one;
		else
			return serialize_duplicate(one);
    }
    check_ptr (is_serialized(one) && is_serialized(two));
   
    // calculate new information length
    h_one = (struct shead *) one;
    size_one = h_one->bytelen + (sizeof(uint32_t) * h_one->tuplelen);
    
    h_two = (struct shead *) two;
    size_two= h_two->bytelen + (sizeof(uint32_t) * h_two->tuplelen);

    // allocate space
    out = myalloc (sizeof(struct shead) + size_one + size_two);
    if (!out)
	return NULL;

    // copy information
    h_out = (struct shead *) out;
    h_out->tuplelen = h_one->tuplelen + h_two->tuplelen;
    h_out->bytelen = h_one->bytelen + h_two->bytelen;
   
    if (!memcpy(&out[sizeof(struct shead)],
	    &one[sizeof(struct shead)],size_one))
    	goto cleanup;
    if (!memcpy(&out[sizeof(struct shead) + size_one],
	    &two[sizeof(struct shead)],size_two))
	goto cleanup;
    
    // destroy old information
    if (del){
	    myfree (one);
	    myfree (two);
    }

    return out;

cleanup:
    myfree (out);
    return NULL;
}

#ifdef __KERNEL__
EXPORT_SYMBOL(serialize);
EXPORT_SYMBOL(serialize_add);
EXPORT_SYMBOL(is_serialized);
EXPORT_SYMBOL(serialized_data);
EXPORT_SYMBOL(serialized_totallen);
EXPORT_SYMBOL(serialize_merge);
#endif

