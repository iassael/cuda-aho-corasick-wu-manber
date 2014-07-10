// string.h
// standard string functionality that is not always available
//
// (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// LGPL license applies

#include "macros.h"

#ifdef __KERNEL__
#include <linux/slab.h>
#include <linux/string.h>
#include <linux/types.h>
#else
#include <string.h>
#include <stdint.h>
#endif

#ifdef __KERNEL__
long 	      strtol(const char *in, char **out, int base);
unsigned long strtoul(const char *in, char **out, int base);
unsigned long strtoul_ex(const char *in, char **out, int base, int *err);
char * 	      strdup(const char *in);
#endif // __KERNEL__
uint32_t      strtohost(const char *string, uint16_t *port);

// the following *should* not be here, but strnlen is sometimes missing 
#ifndef strnlen
#define mystrnlen(a,b) ((strlen(a) > b) ? (b) : strlen(a))
#endif

