/*
 * UpdateModule.cpp
 *
 *  Created on: 15 груд. 2010
 *      Author: sanyi
 */

#include "../../include/update/UpdateModule.h"
#include "../../include/datamodel/DataModel.h"
#include "../../include/datamodel/DataModel_PDM.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/optimization.h"

namespace pssalib
{
	namespace update
	{
		UpdateModule::UpdateModule()
		{
		}

		UpdateModule::~UpdateModule()
		{
		}

		////////////////////////////////
		// Update module methods

		//! Schedule delayed reaction
		bool UpdateModule::schedule_delayed(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();

			// Store the reaction
			pssalib::datamodel::DataModel::DelayedReaction
				reaction(ptrData->mu, ptrSimInfo->dTimeSimulation + ptrData->arSubvolumes[ptrData->sv].arDelay[ptrData->mu]);
			ptrData->vQueuedReactions.push_back(reaction);
			std::sort(ptrData->vQueuedReactions.begin(), ptrData->vQueuedReactions.end());

			return true;
		}

		//! Determine update mode & perform update
		bool UpdateModule::do_update(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{ 
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();
			
			// Set update mode
			if(!ptrSimInfo->getOverrideUpdate()) // is the mode set externally?
			{
				if(ptrData->IsDiffusionReaction(ptrData->mu))
				{
				}
				else
				{
					if(pssalib::datamodel::DataModel::D2 == ptrData->arSubvolumes[ptrData->sv].arReactClass[ptrData->mu])
					{
						schedule_delayed(ptrSimInfo);	// Schedule a D2 reaction
						ptrSimInfo->uUpdateMode = 1;	// Update reactants (D2 delayed reaction)
					}
					else if(pssalib::datamodel::DataModel::D0 == ptrData->arSubvolumes[ptrData->sv].arReactClass[ptrData->mu])
						ptrSimInfo->uUpdateMode = 0;	// Update ALL (no delay)
					else if(pssalib::datamodel::DataModel::D1 == ptrData->arSubvolumes[ptrData->sv].arReactClass[ptrData->mu])
					{
						schedule_delayed(ptrSimInfo);	// Schedule a D1 reaction
						return true;					// Update ALL when the reaction is scheduled (D1 delay reaction)
					}

					else
					{
						ptrSimInfo->reportError("UpdateModule::do_update() - Error: unknown reaction type");
						return false;
					}
				}
			}
			
			// Update species population
			if(!this->update_species(ptrSimInfo))
				return false;

			// Update data structures
			if(!this->update_structures(ptrSimInfo))
				return false;

			if(!this->update_volume_structures(ptrSimInfo))
				return false;

			return true;
		}

		//! Update population
		//! Update mode : 0 - ALL, 1 - reactants & store dalyed reaction, 2 - products
		bool UpdateModule::update_species(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();
			pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[ptrData->sv];
			
			if(ptrData->IsDiffusionReaction(ptrData->mu))
			{
				pssalib::datamodel::DataModel::Subvolume &dstsv = ptrData->arSubvolumes[ptrData->svDst];
				UNSIGNED_INTEGER species = ptrData->mu - ptrData->unReactions + 1;
				sv.aruN[species] -= 1;
				dstsv.aruN[species] += 1;
			}
			else
			{
				// Pointer to the stoichiometry matrix
				MatrixVarLen<pssalib::datamodel::DataModel::ChemicalSpecies>	*arNu	= NULL;

				switch(ptrSimInfo->uUpdateMode)
				{
				// Update ALL
				case 0:
					arNu = &ptrData->arU;
					break;
				// Update reactants
				case 1:
					arNu = &ptrData->arUm;
					break;
				// Update products
				case 2:
					arNu = &ptrData->arUp;
					break;
				default:
					// error
					ptrSimInfo->reportError("UpdateModule_PSSACR::update_structures() - Error: invalid update mode.");
					return false;
					break;
				}

				// Fire the reaction
				UNSIGNED_INTEGER nCols = (*arNu).get_cols(ptrData->mu);
				for(UNSIGNED_INTEGER nu = 0; nu < nCols; nu++)
				{
					pssalib::datamodel::DataModel::ChemicalSpecies	* cs = &(*arNu)(ptrData->mu,nu);
					sv.aruN[cs->index] += cs->coefficient;
					
					if (sv.aruN[cs->index] < 0)
					{
						std::cout << "\n\n NOOOOOOOOOOOO 5 " << cs->coefficient << "     " << sv.aruN[cs->index] << " \n\n";
						pssalib::datamodel::DataModel_PDM* ptrData = (pssalib::datamodel::DataModel_PDM *)ptrSimInfo->getDataModel();
						pssalib::datamodel::DataModel_PDM::Subvolume_PDM &sv = ptrData->arSubvolumes_PDM[ptrData->sv];
						std::cout << sv.arPi << "    " << ptrData->sv << "\n\n\n";
						exit(1);
					}
				}

				// Make sure the reservoir remains constant
				sv.aruN[0] = 1;
			}

			return true;
		}

		bool UpdateModule::update_volume_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();

			if (ptrData->unSubvolumes > 1) 
			{
				// Update source volume propensity.
				pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[ptrData->sv];
			
				UNSIGNED_INTEGER k;
				
				k = floor_log2((UNSIGNED_INTEGER)(sv.dTotalPropensity / ptrData->crsdVolume.minValue)) + 1;
				
				//UNSIGNED_INTEGER k = (UNSIGNED_INTEGER)std::floor(fabs(log2(sv.dTotalPropensity / ptrData->crsdVolume.minValue))) + 1;
						
				ptrData->crsdVolume.updateValue(k, ptrData->sv, sv.dTotalPropensity);

				if (ptrData->IsDiffusionReaction(ptrData->mu))
				{
					// Update destination volume propensity.
					pssalib::datamodel::DataModel::Subvolume &dstsv = ptrData->arSubvolumes[ptrData->svDst];
					//UNSIGNED_INTEGER k = floor_log2((UNSIGNED_INTEGER)(sv.dTotalPropensity / ptrData->crsdVolume.minValue)) + 1;
					UNSIGNED_INTEGER k = (UNSIGNED_INTEGER)std::floor(fabs(log2(dstsv.dTotalPropensity / ptrData->crsdVolume.minValue))) + 1;
					ptrData->crsdVolume.updateValue(k, ptrData->svDst, dstsv.dTotalPropensity);
				}
			}

			return true;
		}
	}	// namespace update
}	// namespace pssalib
