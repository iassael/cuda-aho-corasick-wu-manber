/**
 * support/atomic.[ch]
 * streamline wrapper around atomic operations
 *
 * (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
 * email at wdebruij _-at-_ users DOT sourceforge DOT net
 *
 * Based on the original SUNRPC implementation as found in GLIBC.
 * That version follows an MIT-like license.
 * Here LGPL applies.
 * */

#ifdef __KERNEL__
#include <asm/atomic.h>
#else

/** this is obviously NOT atomic.
 *  TODO: fix. at least now we have the calls in place */

#define atomic_t int

#define ATOMIC_INIT(x) (x)

#define atomic_read(x) (*x)
#define atomic_inc(x) ((*x))++
#define atomic_dec(x) ((*x)--)
#define atomic_inc_and_test(x) ( ++(*x) )

#endif

