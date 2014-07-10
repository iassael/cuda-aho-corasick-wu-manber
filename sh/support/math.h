/* support/math.h 
 * mathematical support routines
 *
 * (c) 2009, Willem de Bruijn, Vrije Universiteit Amsterdam
 * GPLv2 applies
 */

#ifndef STREAMLINE_SUPPORT_MATH_H
#define STREAMLINE_SUPPORT_MATH_H

static inline int 
order_log2(unsigned long in)
{
	unsigned long value = in;
	int i, bytelen, order = 0;

	if (in == 0)
		return 0;

	bytelen = sizeof(unsigned long) * 8;

	for (i = 0; i < bytelen; i++) {
		if (value & 0x1)
			order = i;
		value = value >> 1;
	}
	
	// round up
	if (1 << order == in)
		return order;
	else
		return order + 1;
}

#endif /* STREAMLINE_SUPPORT_MATH_H */

