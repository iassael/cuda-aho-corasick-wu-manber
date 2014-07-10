// stack.h
// very simple stack that used to be part of macros.h
//
// (c) 2008, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_AT_- computer.org
//
// BSD License applies

#ifndef SL_SUPPORT_STACK_H
#define SL_SUPPORT_STACK_H

#include "macros.h"

#define MAGICKEY (-1)

// STACK
//
// NB: MAGICKEY is NOT an acceptable value 
// because stack_pop returns this on an empty stack.
//
// a static stack with private length. can be initialized as 
// full or empty
#define __STACK_EX(keyword, intype, inname, inlen, infull)	\
	keyword int stack_##inname##_len = inlen;		\
	keyword int stack_##inname##_filled = infull;		\
	keyword intype stack_##inname##_entries[inlen];		\
								\
	__attribute__((unused)) 				\
	keyword void 						\
	stack_##inname##_clear(int clearfill) {			\
		bzero(stack_##inname##_entries, sizeof(intype) * stack_##inname##_len);	\
		if (clearfill)								\
			stack_##inname##_filled = 0;					\
	}										\
	keyword inline int 								\
	stack_##inname##_push(intype elem) {						\
		if (likely(stack_##inname##_filled < stack_##inname##_len)) {  		\
			stack_##inname##_entries[stack_##inname##_filled++] = elem; 	\
			return 0;							\
		}									\
		else {									\
			dprintf("stack " #inname " overflow\n"); 			\
			return -1;							\
		}									\
	}										\
					\
	keyword intype			\
	stack_##inname##_pop(void) {	\
		if (stack_##inname##_filled) {						\
			return stack_##inname##_entries[--stack_##inname##_filled];	\
		}									\
		else{									\
			dprintf("stack " #inname " underflow\n");			\
			return (intype) MAGICKEY;					\
		}									\
	}	


#define STATIC_STACK(type, name, len, full) __STACK_EX(static, type, name, len, full)
#define STACK(type, name, len, full) __STACK_EX( , type, name, len, full)

#define stack_clear(inname, fill) stack_##inname##_clear(fill)
#define stack_empty(inname) (unlikely(stack_##inname##_filled == 0))
#define stack_push(inname, elem) stack_##inname##_push(elem)
#define stack_pop(inname) stack_##inname##_pop()

/// It may seem complex, with the stack_empty tests, 
//  but that is only to avoid 'underflow' warnings.
#define stack_foreach(inname, elem) \
	for (elem = (stack_empty(inname) ? ((typeof(elem)) MAGICKEY) : stack_pop(inname));\
	     elem != ((typeof(elem)) MAGICKEY);	\
	     elem = (stack_empty(inname) ? ((typeof(elem)) MAGICKEY) : stack_pop(inname)))

#endif /* SL_SUPPORT_STACK_H */

