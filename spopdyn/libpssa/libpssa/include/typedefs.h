/*
 * typedefs.h
 *
 *  Created on: 2 груд. 2010
 *      Author: sanyi
 */

#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

#include "./stdheaders.h"
#include "./oscompat.h"

///////////////////////////////////////////////////////////////
// Debugging
//#define _DEBUG

#ifdef _DEBUG
	#ifndef _WIN32
		#warning "Debug build"
	#endif
#endif
///////////////////////////////////////////////////////////////

#define PSSA_CR_MAX_ITER 100

///////////////////////////////////////////////////////////////
// Define some basic types
typedef	short				SHORT;
typedef	int					INTEGER;
typedef unsigned int		UNSIGNED_INTEGER;
#ifdef __MACH__
	typedef uint64_t			UINT64;
#endif
typedef	double				REAL;
typedef long double			REAL_EXT;

#if defined(__linux__) || defined(__MACH__)
typedef std::string			STRING;
#elif defined(_WIN32)
	#ifdef _UNICODE
		typedef std::wstring	STRING;
	#else
		typedef std::string		STRING;
	#endif
#endif

///////////////////////////////////////////////////////////////
// Types specific for PSSAlib
namespace pssalib
{
	namespace datamodel
	{
		class DataModel;
	}

	typedef	void	(*FCN_REPORTPROGRESS_CALLBACK)	(UNSIGNED_INTEGER, UNSIGNED_INTEGER, void*);
	typedef	void	(*FCN_REPORTERROR_CALLBACK)		(STRING, void*);
	typedef void	(*FCN_REACTION_CALLBACK)		(pssalib::datamodel::DataModel*, REAL, void*);
}

///////////////////////////////////////////////////////////////
// Global functions specific for PSSAlib
namespace pssalib
{
	REAL				getHmu(UNSIGNED_INTEGER	n, UNSIGNED_INTEGER	m);
	REAL				getCmu(REAL volume, Reaction *reaction);


	UNSIGNED_INTEGER    factorial(UNSIGNED_INTEGER n);
}

///////////////////////////////////////////////////////////////
// MatrixVarLen class - matrix with varying row length

template< typename A >
class MatrixVarLen
{
protected:
	unsigned int uRows, uInc;
	unsigned int *uCols, *uCols_alloc;
	A** data;

public:
	//
	// Default constructor
	MatrixVarLen<A> () : uRows(0), uInc(0), uCols(NULL), uCols_alloc(NULL), data(NULL) {  };

	//
	// Copy constructor
	MatrixVarLen<A> (MatrixVarLen<A> &m) :
			uRows(m.uRows),
			uInc(m.uInc),
			uCols(m.uCols),
			uCols_alloc(m.uCols_alloc),
			data(m.data)
	{
		m.data = NULL;
		m.uRows = 0;
		m.uCols = NULL;
		m.uCols_alloc = NULL;
		m.uInc = 0;
	};

	//
	// Creates a varying row length(C[i]) matrix
	MatrixVarLen<A> (unsigned int uM, unsigned int* uC) : uRows(uM) { resize(uM, uC); };

	//
	// Creates an empty varying row length matrix, preallocating n elements for each row
	MatrixVarLen<A> (unsigned int uM, unsigned int uC) : uRows(uM) { reserve(uM, uC); };

	//
	// Destructor
	~MatrixVarLen<A> ()
	{
		free();
	}

	//
	// Get matrix element
	// i, j - index of the element
	inline const A operator()(unsigned int i, unsigned int j) const 	throw(std::runtime_error)
	{
		if((i < uRows)&&(j < uCols[i]))
			return (data[i])[j];
		else
			throw std::runtime_error("MatrixVarLen<A>::operator() - subscript out of range.");
	};

	//
	// Get matrix element
	// i, j - index of the element
	inline A& operator()(unsigned int i, unsigned int j)				throw(std::runtime_error)
	{
		if((i < uRows)&&(j < uCols[i]))
			return (data[i])[j];
		else
			throw std::runtime_error("MatrixVarLen<A>::operator() - subscript out of range.");
	};

	//
	// Adds a new element at the end of the row
	// i - index of the row, elem - element
	void push_back(unsigned int i, A& elem) throw(std::runtime_error)
	{
		if(i < uRows)
		{
			unsigned int uLen = uCols[i];
			if(uLen < uCols_alloc[i])
			{
				// simply add it
				(data[i])[uLen] = elem;
				uCols[i]++;
			}
			else
			{
				// need to allocate more memory
				A * temp = data[i];

				uCols_alloc[i] += uInc;
				data[i] = new A[uCols_alloc[i]];
				if(uLen > 0)
				{
					memcpy(data[i], temp, uLen*sizeof(A));
					delete [] temp;
				}
				(data[i])[uCols[i]++] = elem;
			}
		}
		else
			throw std::runtime_error("MatrixVarLen<A>::push_back() - subscript out of range.");
	}

	// Resize the matrix
	void resize(unsigned int uM, unsigned int* uC)						throw(std::runtime_error)
	{
		if((0 == uM)||(NULL == uC))
			throw std::runtime_error("MatrixVarLen<A>::resize() - invalid arguments.");

		if(data != NULL)
		{
			for(unsigned short  i = 0; i < uRows; i++)
				delete [] (data[i]);
			delete [] data;
			delete [] uCols;
			delete [] uCols_alloc;
		}

		uRows = uM;
		data = new A*[uM];

		uCols = new unsigned int[uM]; uCols_alloc = new unsigned int[uM];
		memcpy(uCols, uC, uM*sizeof(unsigned int));
		memcpy(uCols_alloc, uCols, uM*sizeof(unsigned int));

		//
		// Calculate the total number of columns
		unsigned int inc = 0;
		for(unsigned int  i = 0; i < uM; i++)
		{
			data[i] = new A[uC[i]]; inc += uC[i];
		}

		// Calculate the growth factor
		uInc = inc / uM;
		if(0 == uInc) uInc = 1;
	};

	// Resize the matrix
	void resize(unsigned int uM, unsigned int uC)						throw(std::runtime_error)
	{
		if((0 == uM)||(0 == uC))
			throw std::runtime_error("MatrixVarLen<A>::resize() - invalid arguments.");

		if(data != NULL)
		{
			for(unsigned short  i = 0; i < uRows; i++)
				delete [] (data[i]);
			delete [] data;
			delete [] uCols;
			delete [] uCols_alloc;
		}

		uRows = uM;
		data = new A*[uM];

		uCols = new unsigned int[uM]; uCols_alloc = new unsigned int[uM];
		for(UNSIGNED_INTEGER i = 0; i < uM; i++)
			uCols[i] = uC;
		memcpy(uCols_alloc, uCols, uM*sizeof(unsigned int));

		//
		// Calculate the total number of columns
		for(unsigned int  i = 0; i < uM; i++)
			data[i] = new A[uC];

		// Calculate the growth factor
		uInc = uC;
	};

	// Preallocate memory (constant row lengths)
	void reserve(unsigned int uM, unsigned int nC)						throw(std::runtime_error)
	{
		if(0 == nC)	nC = 1;
		uInc = nC;

		if(data != NULL)
		{
			// Check if more rows are requested
			if(uRows < uM)
			{
				A **				temp_data = data;
				UNSIGNED_INTEGER *	temp_uCols = uCols, * temp_uCols_alloc = uCols_alloc;

				// Adjust the container
				data = new A*[uM];
				// Adjust & fill helper variables
				uCols = new unsigned int[uM]; uCols_alloc = new unsigned int[uM];

				// Nullify the new arrays
				memset(uCols, (unsigned char)0, uM*sizeof(unsigned int));
				memset(uCols_alloc, (unsigned char)0, uM*sizeof(unsigned int));

				// Copy old values
				memcpy(uCols, temp_uCols, uRows*sizeof(unsigned int));
				memcpy(uCols_alloc, temp_uCols_alloc, uRows*sizeof(unsigned int));

				// Copy the data
				for(unsigned int  i = 0; i < uRows; i++)
				{
					data[i] = temp_data[i];
				}


				uRows = uM;

				// Delete temporary variables
				delete [] temp_data;
				delete [] temp_uCols;
				delete [] temp_uCols_alloc;
			}
			else if(uRows > uM)
				throw std::runtime_error("MatrixVarLen<A>::reserve() - invalid arguments.");

			UNSIGNED_INTEGER	uLen;
			A *					temp;
			for(unsigned int  i = 0; i < uRows; i++)
			{
				uLen = uCols_alloc[i];
				if(uLen < nC)
				{
					temp = data[i];
					data[i] = new A[nC];
					if(uLen > 0)
					{
						memcpy(data[i], temp, uLen*sizeof(A));
						delete [] temp;
					}
					uCols_alloc[i] = nC;
				}
			}
		}
		else
		{
			uRows = uM;
			data = new A*[uM];
			uCols = new unsigned int[uRows]; uCols_alloc = new unsigned int[uRows];
			for(unsigned int  i = 0; i < uRows; i++)
			{
				data[i]			= new A[nC];
				uCols[i]		= 0;
				uCols_alloc[i]	= nC;
			}
		}
	};

	// Clear matrix
	inline void clear()
	{
		for(unsigned  int	i = 0; i < uRows; i++)
		{
			for(unsigned  int	j = 0; j < uCols[i]; j++)
				(&(data[i])[j])->A::~A();

			memset(data[i], (unsigned char)0, uCols[i]*sizeof(A));
		}
		memset(uCols, (unsigned char)0, uRows*sizeof(unsigned int));
	};

	// Free resources
	inline void free()
	{
		if(data != NULL)
		{
			for(unsigned int  i = 0; i < uRows; i++)
				delete [] (data[i]);

			uRows = 0;	uInc = 0;

			delete [] data;
			data = NULL;

			delete [] uCols;
			uCols = NULL;

			delete [] uCols_alloc;
			uCols_alloc = NULL;
		}
	}

	// Get number of rows
	inline unsigned int get_rows() const
	{
		return uRows;
	};

	// Get number of columns
	inline unsigned int get_cols(unsigned int row) const
	{
		return uCols[row];
	};

	// Operator<< for console output
	friend std::ostream & operator<<( std::ostream & output, const MatrixVarLen<A> & M)
	{
		output << "Matrix : dimensions[" << M.uRows << ',' << "*]" << std::endl << std::setprecision(7);
		for(unsigned int i = 0; i < M.uRows; i++)
		{
			for(unsigned int j = 0; j < M.uCols[i]; j++)
			{
				output << std::setw(9) << M(i,j) << ' ';
			}
			output << std::endl;
		}
		return output;
	};
};

///////////////////////////////////////////////////////////////
// Class for storing the result of calculation
class DataSorter
{
protected:
	bool				bSortByIndex;
	UNSIGNED_INTEGER	unIndex;
	std::stringstream	ssTemp;

public:
	DataSorter():
		bSortByIndex(false), unIndex(0)
	{

	}

	DataSorter(const DataSorter& ds):
		bSortByIndex(ds.bSortByIndex), unIndex(ds.unIndex)
	{
		ssTemp.str("");
	}

	void	SetIndex(UNSIGNED_INTEGER unI)	{	bSortByIndex = true; unIndex = unI;	};
	void	UnsetIndex()					{	bSortByIndex = false;	};

	// Converts the input array to a string according current sorting settings
	STRING	getString(std::vector<UNSIGNED_INTEGER> a)
	{
		ssTemp.str("");
		if(bSortByIndex)
			ssTemp << std::setw(12) << a[unIndex] << std::setw(1) << ' ';
		else
		{
			for(UNSIGNED_INTEGER i = 0; i < a.size(); i++)
				ssTemp << std::setw(12) << a[i] << std::setw(1) << ' ';
		}
		return ssTemp.str();
	};

	// Comparison (==)
	bool	equals(std::vector<UNSIGNED_INTEGER> a, std::vector<UNSIGNED_INTEGER> b)
	{
		if(a.size() == b.size())
		{
			if(bSortByIndex)
				return (a[unIndex] == b[unIndex]);
			else
			{
				for(UNSIGNED_INTEGER i = 0; i < a.size(); i++)
				{
					if(a[i] != b[i])
						return false;
				}
				return true;
			}
		}
		else
			return false;
	};
	
	// Comparison operator for std::sort
	bool	operator()(std::vector<UNSIGNED_INTEGER> a, std::vector<UNSIGNED_INTEGER> b)	throw(std::runtime_error)
	{
		if(a.size() == b.size())
		{
			if(bSortByIndex)
			{
				return (a[unIndex] < b[unIndex]);
			}
			else
			{
				UNSIGNED_INTEGER sz = a.size() - 1;
				for(UNSIGNED_INTEGER i = 0; i < sz; i++)
				{
					if(a[i] != b[i])
						return (a[i] < b[i]);
				}
				return (a[sz] < b[sz]);
			}
		}
		else
			throw std::runtime_error("DataSorter::operator<() - incomparable N-tuples: incompatible tuple size.");
	};
};


#endif
