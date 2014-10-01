/*
 * GroupingModule_DM.cpp
 *
 *  Created on: 14 груд. 2010
 *      Author: sanyi
 */

#include "../../include/grouping/GroupingModule_DM.h"
#include "../../include/datamodel/DataModel_DM.h"

namespace pssalib
{
	namespace grouping
	{
		GroupingModule_DM::GroupingModule_DM()
		{
		}

		GroupingModule_DM::GroupingModule_DM(GroupingModule &g) :
				GroupingModule(g)
		{
		}

		GroupingModule_DM::~GroupingModule_DM()
		{
		}

		////////////////////////////////
		// Grouping module methods

		//! Break-up non-elementary reactions?
		bool GroupingModule_DM::do_breakup() const
		{
			return	false;
		}

		//! Initialize data structures (called before each trial)
		bool GroupingModule_DM::initialize(pssalib::datamodel::DataModel* ptrData)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_DM* ptrDMData = static_cast<pssalib::datamodel::DataModel_DM*>(ptrData);

			// Call the baseclass method
			if(!GroupingModule::initialize(ptrData))
				return false;

			// Reinitialize data structures
			if(ptrDMData->arSubvolumes_DM != NULL)
				delete [] ptrDMData->arSubvolumes_DM;
			ptrDMData->arSubvolumes_DM = new pssalib::datamodel::DataModel_DM::Subvolume_DM[ptrData->unSubvolumes];
			
			ptrDMData->dTotalPropensity = 0.0;
			for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; ++si)
			{
				pssalib::datamodel::DataModel::Subvolume &sv = ptrDMData->arSubvolumes[si];
				pssalib::datamodel::DataModel_DM::Subvolume_DM &svdm = ptrDMData->arSubvolumes_DM[si];

				UNSIGNED_INTEGER numpi = ptrData->unReactions + (ptrData->unSpecies - 1);
				svdm.arPi = new REAL[numpi];
				memset(svdm.arPi, 0, sizeof(REAL)*numpi);
				sv.dTotalPropensity = 0.0;

				// Fill the propensities array
				// Normal reactions
				pssalib::datamodel::DataModel::ChemicalSpecies *cs;
				for(UNSIGNED_INTEGER mu = 0, nCols = 0; mu < ptrDMData->unReactions; mu++)
				{
					REAL temp = sv.ardC[mu];
					nCols = ptrData->arUm.get_cols(mu);
					for(UNSIGNED_INTEGER nu = 0; nu < nCols; nu++)
					{
						cs = &ptrDMData->arUm(mu,nu);
						temp *= getHmu(sv.aruN[cs->index], abs(cs->coefficient));
					}
					svdm.arPi[mu] = temp;
					sv.dTotalPropensity += temp;
					ptrData->dTotalPropensity += temp;
				}
				
				// Diffusion reactions
				for(UNSIGNED_INTEGER nu = 0; nu < (ptrData->unSpecies - 1); nu++)
				{
					// Isotropic, homogeneous, normal diffusion.
					REAL lumpedP = sv.aruN[nu+1] * sv.ardD[nu];
					lumpedP *= sv.unDiffusionReactions;
					svdm.arPi[ptrDMData->unReactions + nu] = lumpedP;
					sv.dTotalPropensity += lumpedP;
					ptrData->dTotalPropensity += lumpedP;
				}
			}

			return true;
		}
	}
}
