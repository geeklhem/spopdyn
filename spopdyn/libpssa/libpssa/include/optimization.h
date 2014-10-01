/*
 * oscompat.h
 *
 *  Created on: 2012-01-22
 *      Author: Bram Thijssen
 */

#ifndef OPTIMIZATION_H_
#define OPTIMIZATION_H_

// Based on http://aggregate.org/MAGIC/
inline unsigned int floor_log2(unsigned int x)
{
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
				
	x >>= 1;

	x -= ((x >> 1) & 0x55555555);
	x = (((x >> 2) & 0x33333333) + (x & 0x33333333));
	x = (((x >> 4) + x) & 0x0f0f0f0f);
	x += (x >> 8);
	x += (x >> 16);
	return(x & 0x0000003f);
}

#endif