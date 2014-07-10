
// serialize.[ch]
// pack/unpack a bunch of parameters
//
// (c) 2005, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// BSD license applies

#ifndef WJDB_SERIALIZE_H
#define WJDB_SERIALIZE_H

#ifdef __KERNEL__
#include <linux/types.h>
#include <linux/string.h>
#else
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#endif

/** serialize a bunch of parameters
 *
 * call this function with a list of (int, char*) tuples, whereby
 * int contains the #bytes that should be copied starting at the pointer.
 *
 * for example: serial_string = serialize(2, sizeof(int), &my_int, 10, "0612345678");
 * NB: for \0 terminated strings, don't forget to serialize strlen+1 bytes, instead of strlen
 * 
 *  
 * @param argcount contains the number of tuples
 * @return a newly allocated memory block containing the serialized structure
 *
 * the function takes platform specific int's as input, but generates uint32_t's
 * for its internal datasize headers.
 * */
char * serialize(int argcount, ...);

/// add an element to an existing serialized string
char * serialize_add(char *serialized, int dlen, const void *data);

/// deserialize a string that was previously encoded with serialize(..) 
int deserialize(char *string, ...);

/// return the number of elements that are encoded in the string
static inline unsigned int serialized_tuplelen(const char *serialized)
{
	if (serialized)
		return (unsigned int) ((uint32_t*) serialized)[0];
	else
		return 0;
}

/// return the number of bytes are encoded (i.e., don't include metadata size in this calculation)
static inline unsigned int serialized_bytelen(const char *serialized)
{
	if (serialized)
		return (unsigned int) ((uint32_t*) serialized)[1];
	else
		return 0;
}

/// get the length in bytes of the entire serialized package
int serialized_totallen(const char *serialized);

/// is this string one of our serialized strings? 
int is_serialized(const char *serialized);

/// create a duplicate
char * serialize_duplicate(char * to);

/// merge two strings. the two inputs will be destroyed.
/// @param del set to 1 to delete the original strings
char * serialize_merge(char * one, char * two, int del);

/// return a pointer into the packet string
/// @param dlen may be NULL, otherwise it contains the length of the element on return
/// counting of elements starts at 0
#define ser_data(ser, elemno) serialized_data(ser, elemno, NULL)
const char * serialized_data(const char *serialized, int elemno, int *dlen);

#define serialized_foreach(serialized, i, data, len) \
	for(i=0; \
		i<serialized_tuplelen(serialized) && \
		(data=serialized_data(serialized,i, &len )) != NULL \
	;i++) 

#define serialized_foreach_paired(serialized, i, key, value, keylen, vallen) \
	for(i=0; \
		(i+1<serialized_tuplelen(serialized)) && \
		(key=serialized_data(serialized,i, &keylen )) != NULL && \
		(value=serialized_data(serialized,i+1, &vallen )) != NULL \
	;i+=2) 

// weak typing

#define SER_STR(x)  strlen(x) + 1, x
#define SER_U16(x)  2, (&x)
#define SER_U32(x)  4, (&x)
#define SER_U64(x)  8, (&x)
#define SER_BYTE(x) 1, (&x)
#define SER_INT(x)  sizeof(int), (&x)
#define SER_LONG(x)  sizeof(long), (&x)
#define SER_TYPED(x, type) sizeof(type), (&x)
#define SER_BLOB(len, x) len, x

#define DSER_STR(str, n)  (char*) ser_data(str, n)
#define DSER_U16(str, n)  *(uint16_t*) ser_data(str, n)
#define DSER_U32(str, n)  *(uint32_t*) ser_data(str, n)
#define DSER_BYTE(str, n) *(char*) ser_data(str, n)
#define DSER_INT(str, n)  *(int*) ser_data(str, n)
#define DSER_LONG(str, n)  *(long*) ser_data(str, n)

struct shead {
	uint32_t tuplelen;	///< number of elements in the stream
	uint32_t bytelen;	///< the length in bytes of the payload, i.e. excluding item headers and global header
};

#endif /* WJDB_SERIALIZE_H */

