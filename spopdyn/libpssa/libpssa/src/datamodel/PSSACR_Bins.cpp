/**
 * \file PSSA_Bins.cpp
 * \brief Implements methods for propensity binning for PSRD-CR.
 *
 * \date   01.06.2011
 * \author sanyi
 */

#include "../../include/stdheaders.h"
#include "../../include/typedefs.h"
#include "../../include/datamodel/PSSACR_Bins.h"

namespace pssalib
{
	namespace datamodel
	{
		////////////////////////////////
		// Bin class

		//! Constructor
		PSSACR_Bin::PSSACR_Bin() : dBinSum(0.0), arunBinEl(NULL), unNumBinEl(0), unCapBinEl(0)
		{
			// Do nothing
		}

		//! Copy constructor
		PSSACR_Bin::PSSACR_Bin(const PSSACR_Bin &b) : dBinSum(0.0), arunBinEl(NULL), unNumBinEl(0), unCapBinEl(0)
		{
			if(NULL != b.arunBinEl)
			{
				resize(b.unCapBinEl);
				unNumBinEl = b.unNumBinEl;
				memcpy(b.arunBinEl, arunBinEl, unNumBinEl*sizeof(UNSIGNED_INTEGER));
				dBinSum	= b.dBinSum;
			}
		}

		//! Destructor
		PSSACR_Bin::~PSSACR_Bin()
		{
			// Free memory
			clear();
		}

		//! Resize the bin
		void 	PSSACR_Bin::resize(UNSIGNED_INTEGER unsize)
		{
			clear();
			if(unsize > 0)
			{
				arunBinEl = new UNSIGNED_INTEGER[unsize];
				unCapBinEl = unsize;
			}
		}

		//! Clear the bin
		void 	PSSACR_Bin::clear()
		{
			if(NULL != arunBinEl)
			{
				delete [] arunBinEl;
				unCapBinEl = 0;
				unNumBinEl = 0;
				dBinSum	   = 0.0;
			}
		}

		//! Add an index to the list
		//! \param el UNSIGNED_INTEGER Propensity index that will be added to this bin
		//! \return Returns the index of the newly added element
		UNSIGNED_INTEGER	PSSACR_Bin::push_back(UNSIGNED_INTEGER el)			throw(std::runtime_error)
		{
			if(unNumBinEl < unCapBinEl)
				arunBinEl[unNumBinEl] = el;
			else
				throw std::runtime_error("PSSACR_Bin::push_back() - number of elements exceeds vector capacity.");
			return unNumBinEl++;
		}

		//! Remove the element at the given index
		//! \attention Constant time!
		//! \param idx UNSIGNED_INTEGER Index of the element
		//! \internal \return Returns index of the binned element that was moved
		UNSIGNED_INTEGER	PSSACR_Bin::remove_at(UNSIGNED_INTEGER idx)			throw(std::runtime_error)
		{
			if( (idx + 1) == unNumBinEl)
				--unNumBinEl;
			else if(idx < unNumBinEl)
				arunBinEl[idx] = arunBinEl[--unNumBinEl];
			else
				throw std::runtime_error("PSSACR_Bin::remove_at() - invalid index supplied as argument.");

			return arunBinEl[idx]; // element number to be updated
		}

		//! Retrieve element at the given index
		//! \param idx UNSIGNED_INTEGER Index of the element
		//! \return Returns the stored propensity index
		UNSIGNED_INTEGER	PSSACR_Bin::get_at(UNSIGNED_INTEGER idx)	const	throw(std::runtime_error)
		{
			if(idx < unNumBinEl)
				return arunBinEl[idx];
			else
				throw std::runtime_error("PSSACR_Bin::get_at() - invalid index in arguments.");
		}

		//! Returns number of entries in index list
		//! \return an UNSIGNED_INTEGER representing number of entries in index list
		UNSIGNED_INTEGER	PSSACR_Bin::size()	const
		{
			return unNumBinEl;
		}

		////////////////////////////////
		// Bin storage & mapping

		//! Default constructor
		PSSACR_Bins::PSSACR_Bins()	:	binVals (NULL), unVals(0)
		{
			// Goggle dense_hash_map specific
#ifdef __USE_GOOGLE_HASH_MAP
			mapBins.set_empty_key(0);
#endif
		}

		//! Destructor
		PSSACR_Bins::~PSSACR_Bins()
		{
			clear();
		}

		/**
		 * \brief Updates value in a bin
		 *
		 * \param bin_no_new UNSIGNED_INTEGER New bin number for the given index. Zero denotes a deletion.
		 * \param idx UNSIGNED_INTEGER Item index (e.g., row index in group sum array).
		 * \param val REAL New value to be set for item at index idx.
		 */
		void PSSACR_Bins::updateValue(UNSIGNED_INTEGER bin_no_new, UNSIGNED_INTEGER idx, REAL val)
		{
			if(val > 0.0)
			{
				if(0 == binVals[idx].bin_no)
				{
					// insert
					it = mapBins.find(bin_no_new);
					if(mapBins.end() == it)
					{
						PSSACR_Bin bin;
						insRes = mapBins.insert(MAP_BINS::value_type(bin_no_new, bin));

						it = insRes.first;
						it->second.dBinSum = val;
						it->second.resize(unVals);

						binVals[idx].idx = it->second.push_back(idx);
						binVals[idx].bin_no = bin_no_new;
						binVals[idx].val = val;
					}
					else
					{
						it->second.dBinSum += val;

						binVals[idx].idx = it->second.push_back(idx);
						binVals[idx].bin_no = bin_no_new;
						binVals[idx].val = val;
					}
				}
				else if(bin_no_new == binVals[idx].bin_no)
				{
					// update
					it = mapBins.find(binVals[idx].bin_no);
					if(mapBins.end() == it)
					{
						std::cerr << "Bin does not exist!\nBin = " << binVals[idx].bin_no << " : idx = " << idx << " : val = " << val << std::endl;
					}
					else
					{
						it->second.dBinSum +=  val - binVals[idx].val;

						binVals[idx].val = val;
					}
				}
				else
				{
					// update & move
					it = mapBins.find(binVals[idx].bin_no);
					if(mapBins.end() == it)
					{
						std::cerr << "Bin does not exist!\nBin = " << binVals[idx].bin_no << " : idx = " << idx << " : val = " << val << std::endl;
					}
					else
					{
						// remove at old position & update the index of swapped element
						binVals[it->second.remove_at(binVals[idx].idx)].idx = binVals[idx].idx;

						it->second.dBinSum -= binVals[idx].val;

						it = mapBins.find(bin_no_new);
						if(mapBins.end() == it)
						{
							// insert new bin
							PSSACR_Bin bin;
							insRes = mapBins.insert(MAP_BINS::value_type(bin_no_new, bin));

							it = insRes.first;
							it->second.dBinSum = val;
							it->second.resize(unVals);

							binVals[idx].idx = it->second.push_back(idx);
						}
						else
						{
							it->second.dBinSum +=  val;

							binVals[idx].idx = it->second.push_back(idx);
						}

						binVals[idx].bin_no = bin_no_new;
						binVals[idx].val = val;
					}
				}
			}
			else
			{
				it = mapBins.find(binVals[idx].bin_no);
				if(mapBins.end() != it)
				{
					// delete empty bin entry & update the index of swapped element
					binVals[it->second.remove_at(binVals[idx].idx)].idx = binVals[idx].idx;
					it->second.dBinSum -= binVals[idx].val;

					binVals[idx].bin_no = 0;
					binVals[idx].val = 0.0;
					binVals[idx].idx = -1;
				}
			}
		}
	}
}
