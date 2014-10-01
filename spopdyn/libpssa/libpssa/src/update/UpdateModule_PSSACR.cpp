/*
 * UpdateModule_PSSACR.cpp
 *
 *  Created on: 22 лют. 2011
 *      Author: sanyi
 */

#include "../../include/datamodel/DataModel_PSSACR.h"
#include "../../include/update/UpdateModule_PSSACR.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/optimization.h"

namespace pssalib
{
	namespace update
	{
		UpdateModule_PSSACR::UpdateModule_PSSACR()
		{
			vIdxSigmaChanged.reserve(20);
		}

		UpdateModule_PSSACR::~UpdateModule_PSSACR()
		{
		}

		////////////////////////////////
		// Update module methods

		//! Update structures
		bool UpdateModule_PSSACR::update_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Call the base class method
			if(!UpdateModule_PDM::update_structures(ptrSimInfo))
				return false;

			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_PSSACR *ptrPSRDCRData = static_cast<pssalib::datamodel::DataModel_PSSACR*>(ptrSimInfo->getDataModel());
			pssalib::datamodel::DataModel::Subvolume &sv = ptrPSRDCRData->arSubvolumes[ptrPSRDCRData->sv];
			pssalib::datamodel::DataModel_PDM::Subvolume_PDM &svpdm = ptrPSRDCRData->arSubvolumes_PDM[ptrPSRDCRData->sv];
			pssalib::datamodel::DataModel_PSSACR::Subvolume_PSSACR &svcr = ptrPSRDCRData->arSubvolumes_PSSACR[ptrPSRDCRData->sv];
			
			if(ptrPSRDCRData->IsDiffusionReaction(ptrPSRDCRData->mu))
			{
				pssalib::datamodel::DataModel_PDM::Subvolume_PDM &dstsvpdm = ptrPSRDCRData->arSubvolumes_PDM[ptrPSRDCRData->svDst];
				pssalib::datamodel::DataModel_PSSACR::Subvolume_PSSACR &dstsvcr = ptrPSRDCRData->arSubvolumes_PSSACR[ptrPSRDCRData->svDst];

				UNSIGNED_INTEGER species = ptrPSRDCRData->mu - ptrPSRDCRData->unReactions + 1;
				
				vIdxSigmaChanged.clear();
				update_structures_species(ptrPSRDCRData, ptrPSRDCRData->sv, species);
				REAL invMinSigma = 1.0 / svcr.crsdSigma.minValue;
				std::vector<UNSIGNED_INTEGER>::iterator	it, itEnd = vIdxSigmaChanged.end();
				for(it = vIdxSigmaChanged.begin();it != itEnd; it++)
				{
					UNSIGNED_INTEGER unIdx = (*it);
					INTEGER iTemp1;
					if(svpdm.arSigma[unIdx] > 0.0)
						iTemp1 = (INTEGER)floor(fabs(log2(svpdm.arSigma[unIdx] * invMinSigma))) + 1;
					else
						iTemp1 = 0;
					svcr.crsdSigma.updateValue(iTemp1, unIdx, svpdm.arSigma[unIdx]);
				}
				
				vIdxSigmaChanged.clear();
				update_structures_species(ptrPSRDCRData, ptrPSRDCRData->svDst, species);
				invMinSigma = 1.0 / dstsvcr.crsdSigma.minValue;
				itEnd = vIdxSigmaChanged.end();
				for(it = vIdxSigmaChanged.begin();it != itEnd; it++)
				{
					UNSIGNED_INTEGER unIdx = (*it);
					INTEGER iTemp1;
					if(dstsvpdm.arSigma[unIdx] > 0.0)
						iTemp1 = (INTEGER)floor(fabs(log2(dstsvpdm.arSigma[unIdx] * invMinSigma))) + 1;
					else
						iTemp1 = 0;
					dstsvcr.crsdSigma.updateValue(iTemp1, unIdx, dstsvpdm.arSigma[unIdx]);
				}
			}
			else
			{
				// Pointer to the stoichiometry matrix
				MatrixVarLen<pssalib::datamodel::DataModel::ChemicalSpecies> *arNu = NULL;

				switch(ptrSimInfo->uUpdateMode)
				{
				// Update ALL
				case 0:
					arNu = &ptrPSRDCRData->arU;
					break;
				// Update reactants
				case 1:
					arNu = &ptrPSRDCRData->arUm;
					break;
				// Update products
				case 2:
					arNu = &ptrPSRDCRData->arUp;
					break;
				default:
					// error
					ptrSimInfo->reportError("UpdateModule_PSSACR::update_structures() - Error: invalid update mode.");
					return false;
					break;
				}

				// Clear temporary variables
				vIdxSigmaChanged.clear();

				// Helper variables
				UNSIGNED_INTEGER	unCols = arNu->get_cols(ptrPSRDCRData->mu);
				for(UNSIGNED_INTEGER nu = 0; nu < unCols; nu++)
				{
					// Get the data
					const pssalib::datamodel::DataModel_PDM::ChemicalSpecies species = (*arNu)(ptrPSRDCRData->mu,nu);
					update_structures_species(ptrPSRDCRData, ptrPSRDCRData->sv, species.index);
				}
				
				REAL invMinSigma = 1.0 / svcr.crsdSigma.minValue;
				std::vector<UNSIGNED_INTEGER>::iterator	it, itEnd = vIdxSigmaChanged.end();
				for(it = vIdxSigmaChanged.begin();it != itEnd; it++)
				{
					UNSIGNED_INTEGER unIdx = (*it);
					INTEGER iTemp1;
					if(svpdm.arSigma[unIdx] > 0.0)
						iTemp1 = (INTEGER)floor(fabs(log2(svpdm.arSigma[unIdx] * invMinSigma))) + 1;
					else
						iTemp1 = 0;
					svcr.crsdSigma.updateValue(iTemp1, unIdx, svpdm.arSigma[unIdx]);
				}
			}

			return true;
		}

		void UpdateModule_PSSACR::update_structures_species(pssalib::datamodel::DataModel_PSSACR* ptrPSRDCRData, UNSIGNED_INTEGER svix, UNSIGNED_INTEGER species)
		{
			pssalib::datamodel::DataModel_PDM::Subvolume_PDM &svpdm = ptrPSRDCRData->arSubvolumes_PDM[svix];
			pssalib::datamodel::DataModel_PSSACR::Subvolume_PSSACR &svcr = ptrPSRDCRData->arSubvolumes_PSSACR[svix];

			UNSIGNED_INTEGER U3_rowlen = ptrPSRDCRData->arU3.get_cols(species);

			vIdxSigmaChanged.push_back(species);

			// For each reaction that contains current species
			for(UNSIGNED_INTEGER l = 0; l < U3_rowlen; l++)
			{
				const pssalib::datamodel::DataModel_PDM::PropensityIndex propIdx = ptrPSRDCRData->arU3(species, l);

				vIdxSigmaChanged.push_back(propIdx.i);

				UNSIGNED_INTEGER iTemp1 = floor_log2((UNSIGNED_INTEGER)(svpdm.arPi(propIdx.i,propIdx.j) / svcr.crsdPi[propIdx.i].minValue)) + 1;
				//INTEGER iTemp1 = (INTEGER)floor(fabs(log2(svpdm.arPi(propIdx.i,propIdx.j) / svcr.crsdPi[propIdx.i].minValue))) + 1;
				svcr.crsdPi[propIdx.i].updateValue(iTemp1, propIdx.j, svpdm.arPi(propIdx.i,propIdx.j));
			}
		}
	}	// namespace update
}	// namespace pssalib
