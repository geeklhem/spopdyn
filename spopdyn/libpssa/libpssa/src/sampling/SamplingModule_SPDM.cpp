/*
 * SamplingModule_PDM.cpp
 *
 *  Created on: 14 груд. 2010
 *      Author: sanyi
 */

#include "../../include/stdheaders.h"
#include "../../include/sampling/SamplingModule_SPDM.h"
#include "../../include/datamodel/DataModel_SPDM.h"
#include "../../include/datamodel/SimulationInfo.h"

namespace pssalib
{
	namespace sampling
	{
		SamplingModule_SPDM::SamplingModule_SPDM()
		{
		}

		SamplingModule_SPDM::~SamplingModule_SPDM()
		{
		}

		////////////////////////////////
		// Sampling module methods

		//! Get next reaction index
		bool SamplingModule_SPDM::sample_reaction(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_SPDM* ptrSPDMData = static_cast<pssalib::datamodel::DataModel_SPDM*>(ptrSimInfo->getDataModel());
			pssalib::datamodel::DataModel::Subvolume &sv = ptrSPDMData->arSubvolumes[ptrSPDMData->sv];
			pssalib::datamodel::DataModel_PDM::Subvolume_PDM &svpdm = ptrSPDMData->arSubvolumes_PDM[ptrSPDMData->sv];
			pssalib::datamodel::DataModel_SPDM::Subvolume_SPDM &svspdm = ptrSPDMData->arSubvolumes_SPDM[ptrSPDMData->sv];

			// Temporary variables
			REAL temp1, temp2;
			UNSIGNED_INTEGER i, j, N = ptrSPDMData->unSpecies, tempI, tempJ = 0;

			// Sample reaction
			temp1 = gsl_rng_uniform_pos (ptrRNG) * sv.dTotalPropensity;
			temp2 = 0.0;
			for(i = 0; i < N; i++)
			{
				temp2 += svpdm.arSigma[svspdm.arIdxSigma[i]];
				// TODO : Which sign should I use here
				if(temp1 < temp2)
					break;
			}
			// Store position
			ptrSPDMData->uI = i;
			tempI = svspdm.arIdxSigma[i];

			temp1 = (temp1 - temp2 + svpdm.arSigma[tempI]) / (REAL)sv.aruN[tempI];
			temp2 = 0.0;
			N = svpdm.arPi.get_cols(tempI);
			for(j = 0; j < N; j++)
			{
				tempJ = svspdm.arIdxPi(tempI,j);
				temp2 += svpdm.arPi(tempI,tempJ);
				if(temp1 < temp2)
					break;
			}
			// Store position
			ptrSPDMData->uJ = j;

			// Set the next reaction index
			ptrSPDMData->mu = ptrSPDMData->aruL(tempI,tempJ);

			return true;
		}
	}	// namespace sampling
}	// namespace pssalib
