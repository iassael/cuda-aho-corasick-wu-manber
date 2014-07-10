// endian.h
// detect and cope with varying endianness
//
// (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// LGPL license applies

#include "macros.h"

#define ENDIAN_BIG	0x0
#define ENDIAN_LITTLE	0x1

__attribute__((pure)) static inline int arch_get_endianness(void)
{
#if defined i386
	return ENDIAN_LITTLE;
#elif defined sparc
	return ENDIAN_BIG;
#elif defined ppc || defined powerpc
	return ENDIAN_BIG;
#elif defined armbe
	return ENDIAN_BIG;
#else
	int16_t one = 1;
	char *cp = (char*)&one;
	if ( *cp == 0 )
		return ENDIAN_LITTLE;
	return ENDIAN_BIG;
#endif
}

#if defined i386 || defined x86-64
#define SL_BYTEORDER ENDIAN_LITTLE
#elif defined sparc
#define SL_BYTEORDER ENDIAN_BIG
#elif defined __ARMEB__
#define SL_BYTEORDER ENDIAN_BIG
#elif defined ppc || defined powerpc
#define SL_BYTEORDER ENDIAN_BIG
#else
#warning "cannot predefine endianness"
#endif

/// some archs (sun) have 8byte pointers but 4 byte ints, then *(int*) will fail
/// use this as alternative
#define swap16(A)  ((((uint16_t)(A) & 0xff00) >> 8) | \
                   (((uint16_t)(A) & 0x00ff) << 8))
#define swap32(A)  ((((uint32_t)(A) & 0xff000000) >> 24) | \
                   (((uint32_t)(A) & 0x00ff0000) >> 8)  | \
                   (((uint32_t)(A) & 0x0000ff00) << 8)  | \
                   (((uint32_t)(A) & 0x000000ff) << 24))

// swap on not equal: if endian is unequal to local endianness then swap
static inline uint16_t swap16_ne(uint16_t var, int endian)
{
	if (likely(endian == arch_get_endianness()))
		return var;
	return swap16(var);
}

