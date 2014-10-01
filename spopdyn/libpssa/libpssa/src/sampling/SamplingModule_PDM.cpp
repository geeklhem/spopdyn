/*
 * SamplingModule_PDM.cpp
 *
 *  Created on: 14 груд. 2010
 *      Author: sanyi
 */

#include "../../include/stdheaders.h"
#include "../../include/sampling/SamplingModule_PDM.h"
#include "../../include/datamodel/DataModel_PDM.h"
#include "../../include/datamodel/SimulationInfo.h"

namespace pssalib
{
	namespace sampling
	{
		SamplingModule_PDM::SamplingModule_PDM()
		{
		}

		SamplingModule_PDM::~SamplingModule_PDM()
		{
		}

		////////////////////////////////
		// Sampling module methods

		//! Get next reaction index
		bool SamplingModule_PDM::sample_reaction(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_PDM* ptrPDMData = static_cast<pssalib::datamodel::DataModel_PDM*>(ptrSimInfo->getDataModel());
			pssalib::datamodel::DataModel::Subvolume &sv = ptrPDMData->arSubvolumes[ptrPDMData->sv];
			pssalib::datamodel::DataModel_PDM::Subvolume_PDM &svpdm = ptrPDMData->arSubvolumes_PDM[ptrPDMData->sv];

			// Temporary variables
			REAL temp1, temp2;
			UNSIGNED_INTEGER i, j, N = ptrPDMData->unSpecies;

			// Sample reaction
			temp1 = gsl_rng_uniform_pos (ptrRNG) * sv.dTotalPropensity;
			temp2 = 0.0;
			for(i = 0; i < N; i++)
			{
				temp2 += svpdm.arSigma[i];
				// TODO : Which sign should I use here
				if(temp1 < temp2)
					break;
			}

			temp1 = (temp1 - temp2 + svpdm.arSigma[i]) / (REAL)sv.aruN[i];
			temp2 = 0.0;
			N = svpdm.arPi.get_cols(i);
			for(j = 0; j < N; j++)
			{
				temp2 += svpdm.arPi(i,j);
				if(temp1 < temp2)
					break;
			}

			// Set the next reaction index
			ptrPDMData->mu = ptrPDMData->aruL(i,j);

			return true;
		}
	}	// namespace sampling
}	// namespace psrdlib
