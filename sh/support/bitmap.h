// bitmap.h
// support for per-bit operations
//
// (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// LGPL license applies

// set a bitmap; counting starts at 0
static inline void bitmap_set(char *bitmap, unsigned int element)
{
	bitmap[element >> 3] |= (1 << (element & 0x7));
}

static inline int bitmap_isset(const char *bitmap, unsigned int element)
{
	return (bitmap[element >> 3] & (1 << (element & 0x7))) ? 1 : 0;
}

static inline void bitmap_clear(char *bitmap, unsigned int element)
{
	bitmap[element >> 3] &= ~(1 << (element & 0x7));
}

