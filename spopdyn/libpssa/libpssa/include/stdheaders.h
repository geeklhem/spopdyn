/*
 * stdheaders.h
 *
 *  Created on: 28 лист. 2010
 *      Author: sanyi
 */

#ifndef STDHEADERS_H_
#define STDHEADERS_H_

#include <iostream> 		// Standard I/O streams
#include <fstream> 			// Standard File I/O streams
#include <iomanip>			// Standard header for console output manipulation
#include <stdio.h>  		// Standard I/O
#include <string.h>			// Standard String header
#include <sstream>			// Standard String Stream Header
#include <stdlib.h> 		// Standard C library

// time related routines
#if defined(__linux__)
	#include <time.h>
#elif defined(__MACH__)
	#include<mach/mach.h>
	#include<mach/mach_time.h>
#elif defined(_WIN32)
	#define NOMINMAX
	#include <windows.h>
#endif

// Standard C++ headers
#include <cassert>			// asserts
#include <cerrno>
#include <stdexcept>

// definition of mkdir
#if defined(__linux__) || defined(__MACH__)
	#include <sys/stat.h>
	#include <sys/types.h>
#elif defined(_WIN32)
	#include <direct.h>
#endif

// libM
#define _USE_MATH_DEFINES	// Need it for M_PI to be defined in cmath.h
#include <cmath>			// basic math library
#include <limits>			// Need it for min & lowest

// GNU Scientific Library
#include <gsl/gsl_rng.h>	// GSL random generators
#include <gsl/gsl_randist.h>

// STL
#include <algorithm>		// STL algorithms
#include <vector>   		// STL vector container

// libSBML
#include <sbml/SBMLTypes.h>	// Systems Biology Markup Language library
LIBSBML_CPP_NAMESPACE_USE

// Define hashmap type for PSSACR_Bins class
//#define __USE_GOOGLE_HASH_MAP

#ifdef __USE_GOOGLE_HASH_MAP
#include <google/dense_hash_map>
#define HASHMAP google::dense_hash_map
#else
#define HASHMAP boost::unordered_map
#endif

// Boost headers
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/log1p.hpp>

#endif /* STDHEADERS_H_ */
