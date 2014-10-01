/*
 * UpdateModule_PDM.cpp
 *
 *  Created on: 15 груд. 2010
 *      Author: sanyi
 */

#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/update/UpdateModule_SPDM.h"
#include "../../include/datamodel/DataModel_SPDM.h"

namespace pssalib
{
	namespace update
	{
		UpdateModule_SPDM::UpdateModule_SPDM()
		{
		}

		UpdateModule_SPDM::~UpdateModule_SPDM()
		{
		}

		////////////////////////////////
		// Update module methods

		//! Update structures
		bool UpdateModule_SPDM::update_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Call the base class method
			if(!UpdateModule_PDM::update_structures(ptrSimInfo))
				return false;

			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_SPDM* ptrSPDMData = static_cast<pssalib::datamodel::DataModel_SPDM*>(ptrSimInfo->getDataModel());
			
			if(ptrSPDMData->IsDiffusionReaction(ptrSPDMData->mu))
			{
			}
			else
			{
				pssalib::datamodel::DataModel::Subvolume &sv = ptrSPDMData->arSubvolumes[ptrSPDMData->sv];
				pssalib::datamodel::DataModel_PDM::Subvolume_PDM &svpdm = ptrSPDMData->arSubvolumes_PDM[ptrSPDMData->sv];
				pssalib::datamodel::DataModel_SPDM::Subvolume_SPDM &svspdm = ptrSPDMData->arSubvolumes_SPDM[ptrSPDMData->sv];

				// Temporary variables
				UNSIGNED_INTEGER temp, tempI = svspdm.arIdxSigma[ptrSPDMData->uI];

				if(0 != ptrSPDMData->uI)
				{
					// Swap with preceding
					temp = svspdm.arIdxSigma[ptrSPDMData->uI];
					svspdm.arIdxSigma[ptrSPDMData->uI] = svspdm.arIdxSigma[ptrSPDMData->uI - 1];
					svspdm.arIdxSigma[ptrSPDMData->uI - 1] = temp;
				}

				if(0 != ptrSPDMData->uJ)
				{
					// Swap with preceding
					temp = svspdm.arIdxPi(tempI,ptrSPDMData->uJ);
					svspdm.arIdxPi(tempI,ptrSPDMData->uJ) = svspdm.arIdxPi(tempI,ptrSPDMData->uJ - 1);
					svspdm.arIdxPi(tempI,ptrSPDMData->uJ - 1) = temp;
				}
			}

			return true;
		}
	}	// namespace update
}	// namespace pssalib
