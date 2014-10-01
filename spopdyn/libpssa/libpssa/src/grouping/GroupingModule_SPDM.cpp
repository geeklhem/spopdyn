/*
 * GroupingModule_PDM.cpp
 *
 *  Created on: 14 груд. 2010
 *      Author: sanyi
 */

#include "../../include/grouping/GroupingModule_SPDM.h"
#include "../../include/datamodel/DataModel_SPDM.h"

namespace pssalib
{
	namespace grouping
	{
		GroupingModule_SPDM::GroupingModule_SPDM()
		{
		}

		GroupingModule_SPDM::GroupingModule_SPDM(GroupingModule &g) :
				GroupingModule_PDM(g)
		{
		}

		GroupingModule_SPDM::~GroupingModule_SPDM()
		{
		}

		////////////////////////////////
		// Grouping module methods

		//! Initialize data structures (called before each trial)
		bool GroupingModule_SPDM::initialize(pssalib::datamodel::DataModel* ptrData)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_SPDM* ptrSPDMData = static_cast<pssalib::datamodel::DataModel_SPDM*>(ptrData);

			// Call the baseclass method
			if(!GroupingModule_PDM::initialize(ptrData))
				return false;

			// Clean the data members
			if(ptrSPDMData->arSubvolumes_SPDM != NULL)
				delete [] ptrSPDMData->arSubvolumes_SPDM;
			ptrSPDMData->arSubvolumes_SPDM = new pssalib::datamodel::DataModel_SPDM::Subvolume_SPDM[ptrData->unSubvolumes];

			for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; ++si)
			{
				pssalib::datamodel::DataModel_PDM::Subvolume_PDM &svpdm = ptrSPDMData->arSubvolumes_PDM[si];
				pssalib::datamodel::DataModel_SPDM::Subvolume_SPDM &svspdm = ptrSPDMData->arSubvolumes_SPDM[si];

				svspdm.arIdxSigma = new UNSIGNED_INTEGER[ptrSPDMData->unSpecies];
				for(UNSIGNED_INTEGER i = 0; i < ptrSPDMData->unSpecies; ++i)
					svspdm.arIdxSigma[i] = svpdm.arPi.get_cols(i);

				svspdm.arIdxPi.resize(ptrSPDMData->unSpecies, svspdm.arIdxSigma);

				for(UNSIGNED_INTEGER i = 0, n; i < ptrSPDMData->unSpecies; ++i)
				{
					n = svspdm.arIdxSigma[i];
					for(UNSIGNED_INTEGER j = 0; j < n; j++)
						svspdm.arIdxPi(i,j) = j;
					svspdm.arIdxSigma[i] = i;
				}
			}

			return true;
		}
	}
}
