/*
 * oscompat.h
 *
 *  Created on: 2011-10-28
 *      Author: Bram Thijssen
 */

#ifndef OSCOMPAT_H_
#define OSCOMPAT_H_

#ifdef _WIN32
	inline double __cdecl log2(double _X) { return log(_X) * 1.442695040888963 /*= 1 / log(2.)*/; }
#endif

#endif
