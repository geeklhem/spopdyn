/*
 * UpdateModule_DM.cpp
 *
 *  Created on: 15 груд. 2010
 *      Author: sanyi
 */

//namespace pssalib
//{
//	// Forward declarations
////	class PSSA;
//
//	namespace datamodel
//	{
//		class SimulationInfo;
//		class DataModel;
//	}
//}

#include "../../include/update/UpdateModule_DM.h"
#include "../../include/datamodel/SimulationInfo.h"

namespace pssalib
{
	// Forward declarations
	class PSSA;

	namespace datamodel
	{
		class SimulationInfo;
	}

	namespace update
	{
		UpdateModule_DM::UpdateModule_DM()
		{
		}

		UpdateModule_DM::~UpdateModule_DM()
		{
		}

		////////////////////////////////
		// Update module methods

		void UpdateModule_DM::update_structures_species(pssalib::datamodel::DataModel_DM* ptrDMData, UNSIGNED_INTEGER svix, UNSIGNED_INTEGER species, INTEGER species_coeff, REAL& dlt_a_0)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel::Subvolume &sv = ptrDMData->arSubvolumes[svix];
			pssalib::datamodel::DataModel_DM::Subvolume_DM &svdm = ptrDMData->arSubvolumes_DM[svix];
		  
			// Clean-up
			
			double old_tot_prop = sv.dTotalPropensity;
			sv.dTotalPropensity = 0.0;
			UNSIGNED_INTEGER numpi = ptrDMData->unReactions + ptrDMData->unSpecies - 1;
			memset(svdm.arPi, 0, sizeof(REAL)*numpi);

			// Fill the propensities array
			// Normal reactions
			REAL											temp;
			pssalib::datamodel::DataModel::ChemicalSpecies	*cs;
			for(UNSIGNED_INTEGER mu = 0, nCols = 0; mu < ptrDMData->unReactions; mu++)
			{
				temp = sv.ardC[mu];
				nCols = ptrDMData->arUm.get_cols(mu);
				for(UNSIGNED_INTEGER nu = 0; nu < nCols; nu++)
				{
					// BT_TODO why is this different from the code in GroupingModule_DM?
					cs = &ptrDMData->arUm(mu,nu);
					temp *= getHmu(sv.aruN[cs->index], abs(cs->coefficient));
					
/*					temp *= sv.aruN[cs->index];
					if(cs->coefficient < -1)
						temp *= ((REAL)(sv.aruN[cs->index] - 1))*0.5;			// elementary reactions only*/
				}
				svdm.arPi[mu] = temp;
				sv.dTotalPropensity += temp;
			}
				
			// Diffusion reactions
			for(UNSIGNED_INTEGER nu = 0; nu < ptrDMData->unSpecies - 1; nu++)
			{
				// Isotropic, homogeneous, normal diffusion.
				REAL lumpedP = sv.aruN[nu+1] * sv.ardD[nu];
				lumpedP *= sv.unDiffusionReactions;
				svdm.arPi[ptrDMData->unReactions + nu] = lumpedP;
				sv.dTotalPropensity += lumpedP;
			}
			
			dlt_a_0 = sv.dTotalPropensity - old_tot_prop;
		}
		
		//! Update structures
		bool UpdateModule_DM::update_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_DM* ptrDMData = static_cast<pssalib::datamodel::DataModel_DM*>(ptrSimInfo->getDataModel());
			pssalib::datamodel::DataModel::Subvolume &sv = ptrDMData->arSubvolumes[ptrDMData->sv];

			if (ptrDMData->IsDiffusionReaction(ptrDMData->mu))
			{
//				pssalib::datamodel::DataModel::Subvolume &sv = ptrPDMData->arSubvolumes[ptrPDMData->sv];
//				pssalib::datamodel::DataModel::Subvolume &dstsv = ptrPDMData->arSubvolumes[ptrPDMData->svDst];
				
				UNSIGNED_INTEGER species = 0/*ptrPDMData->mu - ptrPDMData->unReactions + 1*/;
				
				REAL d1 = 0.0;
				update_structures_species(ptrDMData, ptrDMData->sv, species, -1, d1);
				ptrDMData->dTotalPropensity += d1;
//				sv.dTotalPropensity += d1;

				REAL d2 = 0.0;
				update_structures_species(ptrDMData, ptrDMData->svDst, species, 1, d2);
//				dstsv.dTotalPropensity += d2;
				ptrDMData->dTotalPropensity += d2;

//				ptrDMData->dTotalPropensity += d1 + d2;
			}
			else
			{
				// Pointer to the stoichiometry matrix
				MatrixVarLen<pssalib::datamodel::DataModel::ChemicalSpecies> *arNu = NULL;
				MatrixVarLen<pssalib::datamodel::DataModel::ChemicalSpecies> *arNu_r = NULL;

				switch(ptrSimInfo->uUpdateMode)
				{
				// Update ALL
				case 0:
					arNu = &ptrDMData->arU;
					arNu_r = &ptrDMData->arU_r;
					break;
				// Update reactants
				case 1:
					arNu = &ptrDMData->arUm;
					arNu_r = &ptrDMData->arUm_r;
					break;
				// Update products
				case 2:
					arNu = &ptrDMData->arUp;
					arNu_r = &ptrDMData->arUp_r;
					break;
				default:
					// error
					ptrSimInfo->reportError("UpdateModule_PDM::update_structures() - Error: invalid update mode.");
					return false;
					break;
				}
				
				// Check if is not an elementar reaction
				
				REAL dlt_a_0 = 0.0;
				for(UNSIGNED_INTEGER nu = 0; nu < arNu_r->get_cols(ptrDMData->mu); nu++)
				{
					// Get the data
					const pssalib::datamodel::DataModel_DM::ChemicalSpecies species = (*arNu_r)(ptrDMData->mu,nu);

					// Skip reservoir
					if(species.index == 0)
						continue;

					update_structures_species(ptrDMData, ptrDMData->sv, species.index, species.coefficient, dlt_a_0);
					ptrDMData->dTotalPropensity += dlt_a_0;
				}

				// Update total propensity
//				ptrPDMData->dTotalPropensity += dlt_a_0;
//				sv.dTotalPropensity += dlt_a_0;
			}

			return true;
		}
	}	// namespace update
}	// namespace pssalib
