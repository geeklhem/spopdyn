/*
 * SamplingModule_PSSACR.cpp
 *
 *  Created on: 22 лют. 2011
 *      Author: sanyi
 */

#include "../../include/stdheaders.h"
#include "../../include/sampling/SamplingModule_PSSACR.h"
#include "../../include/datamodel/DataModel_PSSACR.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/datamodel/PSSACR_Bins.h"

namespace pssalib
{
	namespace sampling
	{
		SamplingModule_PSSACR::SamplingModule_PSSACR()
		{
		}

		SamplingModule_PSSACR::~SamplingModule_PSSACR()
		{
		}

		////////////////////////////////
		// Sampling module methods

		//! Get next sample
		bool SamplingModule_PSSACR::sample_reaction(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_PSSACR* ptrPSRDCRData = static_cast<pssalib::datamodel::DataModel_PSSACR*>(ptrSimInfo->getDataModel());
			pssalib::datamodel::DataModel::Subvolume &sv = ptrPSRDCRData->arSubvolumes[ptrPSRDCRData->sv];
			pssalib::datamodel::DataModel_PDM::Subvolume_PDM &svpdm = ptrPSRDCRData->arSubvolumes_PDM[ptrPSRDCRData->sv];
			pssalib::datamodel::DataModel_PSSACR::Subvolume_PSSACR &svcr = ptrPSRDCRData->arSubvolumes_PSSACR[ptrPSRDCRData->sv];

			REAL r = 0.0;
			UNSIGNED_INTEGER sI = 0;
			bool success = crSampler.Sample(&svcr.crsdSigma, ptrRNG, sv.dTotalPropensity, sI, r);

			if (success) {
				UNSIGNED_INTEGER sJ = 0;
				success = crSampler.Sample(&svcr.crsdPi[sI], ptrRNG, svpdm.arLambda[sI], sJ, r);

				if (success) {
					ptrPSRDCRData->mu = ptrPSRDCRData->aruL(sI, sJ);
				}
			}

			if (!success) {
				ptrSimInfo->reportError("SamplingModule_PSSACR::sample_reaction - Error: sampling did not converge in given number of iterations.");
			}
			return success;
		}
	}	// namespace sampling
}	// namespace pssalib
