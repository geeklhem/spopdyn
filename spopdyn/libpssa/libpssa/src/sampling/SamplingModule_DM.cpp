/*
 * SamplingModule_DM.cpp
 *
 *  Created on: 14 груд. 2010
 *      Author: sanyi
 */

#include "../../include/stdheaders.h"
#include "../../include/sampling/SamplingModule_DM.h"
#include "../../include/datamodel/DataModel_DM.h"
#include "../../include/datamodel/SimulationInfo.h"

namespace pssalib
{
	namespace sampling
	{
		SamplingModule_DM::SamplingModule_DM()
		{
		}

		SamplingModule_DM::~SamplingModule_DM()
		{
		}

		////////////////////////////////
		// Sampling module methods

		//! Sample reaction index
		bool SamplingModule_DM::sample_reaction(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_DM* ptrDMData = static_cast<pssalib::datamodel::DataModel_DM*>(ptrSimInfo->getDataModel());
			pssalib::datamodel::DataModel::Subvolume &sv = ptrDMData->arSubvolumes[ptrDMData->sv];
			pssalib::datamodel::DataModel_DM::Subvolume_DM &svdm = ptrDMData->arSubvolumes_DM[ptrDMData->sv];

			UNSIGNED_INTEGER mu;
			UNSIGNED_INTEGER M = ptrDMData->unReactions;
			UNSIGNED_INTEGER N = ptrDMData->unSpecies - 1;

			// Sample reaction
			REAL temp1 = gsl_rng_uniform_pos (ptrRNG) * sv.dTotalPropensity;
			REAL temp2 = 0.0;
			for(mu = 0; mu < M + N; mu++)
			{
				temp2 += svdm.arPi[mu];
				if(temp1 <= temp2)
					break;
			}

			ptrDMData->mu = mu;

			return true;
		}
	}	// namespace sampling
}	// namespace psrdlib
