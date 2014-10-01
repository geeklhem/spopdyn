/*
 * GroupingModule_PSSACR.cpp
 *
 *  Created on: 22 лют. 2011
 *      Author: sanyi
 */

#include "../../include/grouping/GroupingModule_PSSACR.h"
#include "../../include/datamodel/DataModel_PSSACR.h"
#include "../../include/datamodel/PSSACR_Bins.h"
#include "../../include/optimization.h"

namespace pssalib
{
	namespace grouping
	{
		GroupingModule_PSSACR::GroupingModule_PSSACR()
		{
		}

		GroupingModule_PSSACR::GroupingModule_PSSACR(GroupingModule &g) :
				GroupingModule_PDM(g)
		{
		}

		GroupingModule_PSSACR::~GroupingModule_PSSACR()
		{
		}

		////////////////////////////////
		// Grouping module methods

		//! Initialize data structures (called before each trial)
		bool GroupingModule_PSSACR::initialize(pssalib::datamodel::DataModel* ptrData)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_PSSACR* ptrPSRDCRData = static_cast<pssalib::datamodel::DataModel_PSSACR*>(ptrData);

			// Temporary variables
			UNSIGNED_INTEGER	k, N;
			REAL				temp1 = 0.0,
								temp2 = 0.0;

			// Call the base class method
			if(!GroupingModule_PDM::initialize(ptrData))
				return false;
			
			if(ptrPSRDCRData->arSubvolumes_PSSACR != NULL)
				delete [] ptrPSRDCRData->arSubvolumes_PSSACR;
			ptrPSRDCRData->arSubvolumes_PSSACR = new pssalib::datamodel::DataModel_PSSACR::Subvolume_PSSACR[ptrData->unSubvolumes];

			for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; ++si)
			{
				////////////////////////////////////////////////
				// Calculate the minimum values
				ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdSigma.minValue = std::numeric_limits<REAL>::max();
				ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi.clear();
				ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi.resize(ptrPSRDCRData->unSpecies);

				for(UNSIGNED_INTEGER i = 0, nPi; i < ptrPSRDCRData->unSpecies; i++)
				{
					nPi = ptrPSRDCRData->arSubvolumes_PDM[si].arPi.get_cols(i);

					temp1 = 0.0;
					if(0 == nPi)
						ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi[i].minValue = 0.0;
					else
						ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi[i].minValue = std::numeric_limits<REAL>::max();

					for(UNSIGNED_INTEGER j = 0; j < nPi; j++)
					{
						temp2 = ptrPSRDCRData->arSubvolumes[si].ardC[ptrPSRDCRData->aruL(i,j)];
						bool self_dep = false;

						for(UNSIGNED_INTEGER l = 0, nU3 = ptrPSRDCRData->arU3.get_cols(i); l < nU3; l++)
						{
							if((ptrPSRDCRData->arU3(i,l).i == i)&&(ptrPSRDCRData->arU3(i,l).j == j))
							{
								temp2 *= 0.5;
								self_dep = true;
								break;
							}
						}
						if(ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi[i].minValue > temp2)
							ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi[i].minValue = temp2;

						if(self_dep)
							temp1 = 2.0 * temp2;
						else
							temp1 = temp2;

						if((temp1 > 0.0)&&(ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdSigma.minValue > temp1))
							ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdSigma.minValue = temp1;
					}
				}

				////////////////////////////////////////////////
				// Compute distribution

				// Initialize the data structures
				ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdSigma.bins.clear();
				ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdSigma.bins.resize(ptrPSRDCRData->unSpecies);

				// Distribute the values across the bins
				for(UNSIGNED_INTEGER i = 0; i < ptrPSRDCRData->unSpecies; i++)
				{
					// Sigma's
					if(0 != ptrPSRDCRData->arSubvolumes_PDM[si].arSigma[i])
					{
						if (ptrData->crsdVolume.minValue == 0.0)
							k = floor_log2((UNSIGNED_INTEGER)(ptrPSRDCRData->arSubvolumes_PDM[si].arSigma[i])) + 1;
						else
							k = floor_log2((UNSIGNED_INTEGER)(ptrPSRDCRData->arSubvolumes_PDM[si].arSigma[i] / ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdSigma.minValue)) + 1;
					  
						//k = (UNSIGNED_INTEGER)floor(abs(log2(ptrPSRDCRData->arSubvolumes_PDM[si].arSigma[i] / ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdSigma.minValue))) + 1;
						ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdSigma.bins.updateValue(k, i, ptrPSRDCRData->arSubvolumes_PDM[si].arSigma[i]);
					}

					// Pi's
					N = ptrPSRDCRData->arSubvolumes_PDM[si].arPi.get_cols(i);
					ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi[i].bins.resize(N);
					for(UNSIGNED_INTEGER j = 0; j < N; j++)
					{
						if(0 == ptrPSRDCRData->arSubvolumes_PDM[si].arPi(i,j))
							continue;
						
						k = floor_log2((UNSIGNED_INTEGER)(ptrPSRDCRData->arSubvolumes_PDM[si].arPi(i,j) / ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi[i].minValue)) + 1;
						//k = (UNSIGNED_INTEGER)std::floor(log2(ptrPSRDCRData->arSubvolumes_PDM[si].arPi(i,j) / ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi[i].minValue)) + 1;
						ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi[i].bins.updateValue(k, j, ptrPSRDCRData->arSubvolumes_PDM[si].arPi(i,j));
					}
				}
			}

			return true;
		}

	}
}
