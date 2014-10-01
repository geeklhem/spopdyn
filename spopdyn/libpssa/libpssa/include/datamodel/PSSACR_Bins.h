/**
 * \file PSSACR_Bins.h
 * \brief Declares methods for propensity binning for PSRD-CR.
 *
 * \details Binning is created along with an instance of \c PSSACR_Bin class.
 * \c PSSACR_Bin class implements the \c updateValue function
 * that handles all the updates necessary when a value of
 * a binned quantity got changed. Check out the source for details.
 *
 * \date   01.06.2011
 * \author sanyi
 */

#ifndef PSSA_BINS_H_
#define PSSA_BINS_H_

#include "../stdheaders.h"
#include "../typedefs.h"

namespace pssalib
{
	namespace datamodel
	{
		////////////////////////////////
		//! \class PSSACR_Bin
		//! \brief Bin class
		class PSSACR_Bin
		{
		////////////////////////////////
		// Attributes
		public:
			//! Sum of bin elements
			REAL				dBinSum;
			//! Indices of binned elements in the original array
			UNSIGNED_INTEGER	*arunBinEl;

			//! \internal auxiliary variables:
			UNSIGNED_INTEGER	unNumBinEl;
			//! \internal auxiliary variables:
			UNSIGNED_INTEGER	unCapBinEl;


		////////////////////////////////
		// Constructors
		public:
			// Constructor
			PSSACR_Bin();

			// Copy constructor
			PSSACR_Bin(const PSSACR_Bin &b);

			// Destructor
	virtual ~PSSACR_Bin();

		////////////////////////////////
		// Methods
		public:
			// Add an index to the list
			UNSIGNED_INTEGER	push_back(UNSIGNED_INTEGER el)			throw(std::runtime_error);
			// Remove the element at the given index
			UNSIGNED_INTEGER	remove_at(UNSIGNED_INTEGER idx)			throw(std::runtime_error);
			// Retrieve element at the given index
			UNSIGNED_INTEGER	get_at(UNSIGNED_INTEGER idx)	const	throw(std::runtime_error);

			// Returns number of entries in index list
			UNSIGNED_INTEGER	size()							const;
			// Resize the bin
			void 				resize(UNSIGNED_INTEGER n);
			// Clear the bin
			void 				clear();
		};

		////////////////////////////////
		//! \class PSSACR_Bins
		//! \brief Bins storage & mapping
		class PSSACR_Bins
		{
		////////////////////////////////
		// Typedefs
		public:
			typedef HASHMAP			<UNSIGNED_INTEGER,	PSSACR_Bin		>					MAP_BINS;
			typedef HASHMAP			<UNSIGNED_INTEGER,	PSSACR_Bin		>::iterator			BINS_ITER;
			typedef HASHMAP			<UNSIGNED_INTEGER,	PSSACR_Bin		>::const_iterator	CONST_BINS_ITER;
			typedef	std::pair		<BINS_ITER,			BINS_ITER		>					PAIR_BINS_ITER;
			typedef	std::pair		<CONST_BINS_ITER,	CONST_BINS_ITER	>					CONST_PAIR_BINS_ITER;

			//! Structure for fast access to bins
			typedef struct tagBinsVals
			{
				//! A copy of binned element value
				REAL						val;
				//! mapping variable: Bin number
				UNSIGNED_INTEGER			bin_no;
				//! mapping variable: Index of the corresponding element in the bin
				UNSIGNED_INTEGER			idx;
			} BinVals;

		////////////////////////////////
		// Attributes
		protected:
			MAP_BINS					mapBins;
			BinVals						*binVals;
			UNSIGNED_INTEGER			unVals;

			//! \internal Temporary variables
			BINS_ITER					it;
			//! \internal Temporary variables
			std::pair<BINS_ITER, bool>	insRes;

		////////////////////////////////
		// Constructors
		public:
			// Constructor
			PSSACR_Bins();

			// Destructor
	virtual ~PSSACR_Bins();

		////////////////////////////////
		// Methods
		public:

			//! Get a pair of iterators
			//! for iterating across all bins
	inline	void getBins(CONST_PAIR_BINS_ITER &pit) const
			{
				pit.first = mapBins.begin();
				pit.second = mapBins.end();
			}

			//! Clear bins
	inline	void clear()
			{
				if(mapBins.size() > 0)
					mapBins.clear();
				if(binVals != NULL)
				{
					delete	[]	binVals;
					binVals = NULL;
					unVals = 0;
				}
			}

			//! Resizes bins, ensuring they have
			//! enough capacity to store <i>N</i> elements
	inline	void resize(UNSIGNED_INTEGER N)
			{
				clear();
				mapBins.rehash(N);
				binVals = new BinVals[N];
				memset(binVals, 0, sizeof(BinVals)*N);
				unVals = N;
			}

			// Updates value in a bin
			void updateValue(UNSIGNED_INTEGER bin_no_new, UNSIGNED_INTEGER idx, REAL val);

	inline	REAL getValue(UNSIGNED_INTEGER idx) const
			{
				if(idx < unVals)
					return binVals[idx].val;
				else
					throw std::runtime_error("PSSACR_Bins::getValue() - invalid index in arguments.");
			}
		};
	}
}

#endif /* PSSA_BINS_H_ */
