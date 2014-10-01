/*
 * UpdateModule_PDM.cpp
 *
 *  Created on: 15 груд. 2010
 *      Author: sanyi
 */

#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/update/UpdateModule_PDM.h"
#include "../../include/datamodel/DataModel_PDM.h"

namespace pssalib
{
	namespace update
	{
		UpdateModule_PDM::UpdateModule_PDM()
		{
		}

		UpdateModule_PDM::~UpdateModule_PDM()
		{
		}

		////////////////////////////////
		// Update module methods

		//! Update structures
		bool UpdateModule_PDM::update_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_PDM* ptrPDMData = static_cast<pssalib::datamodel::DataModel_PDM*>(ptrSimInfo->getDataModel());
			pssalib::datamodel::DataModel::Subvolume &sv = ptrPDMData->arSubvolumes[ptrPDMData->sv];
			pssalib::datamodel::DataModel_PDM::Subvolume_PDM &svpdm = ptrPDMData->arSubvolumes_PDM[ptrPDMData->sv];

			if (ptrPDMData->IsDiffusionReaction(ptrPDMData->mu))
			{
				pssalib::datamodel::DataModel::Subvolume &sv = ptrPDMData->arSubvolumes[ptrPDMData->sv];
				pssalib::datamodel::DataModel::Subvolume &dstsv = ptrPDMData->arSubvolumes[ptrPDMData->svDst];
				
				UNSIGNED_INTEGER species = ptrPDMData->mu - ptrPDMData->unReactions + 1;
				
				REAL d1 = 0.0;
				update_structures_species(ptrPDMData, ptrPDMData->sv, species, -1, d1);
				sv.dTotalPropensity += d1;

				REAL d2 = 0.0;
				update_structures_species(ptrPDMData, ptrPDMData->svDst, species, 1, d2);
				dstsv.dTotalPropensity += d2;

				ptrPDMData->dTotalPropensity += d1 + d2;
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
					arNu = &ptrPDMData->arU;
					arNu_r = &ptrPDMData->arU_r;
					break;
				// Update reactants
				case 1:
					arNu = &ptrPDMData->arUm;
					arNu_r = &ptrPDMData->arUm_r;
					break;
				// Update products
				case 2:
					arNu = &ptrPDMData->arUp;
					arNu_r = &ptrPDMData->arUp_r;
					break;
				default:
					// error
					ptrSimInfo->reportError("UpdateModule_PDM::update_structures() - Error: invalid update mode.");
					return false;
					break;
				}
				
				// Check if is not an elementar reaction
				
				REAL dlt_a_0 = 0.0;
				for(UNSIGNED_INTEGER nu = 0; nu < arNu_r->get_cols(ptrPDMData->mu); nu++)
				{
					// Get the data
					const pssalib::datamodel::DataModel_PDM::ChemicalSpecies species = (*arNu_r)(ptrPDMData->mu,nu);

					// Skip reservoir
					if(species.index == 0)
						continue;

					update_structures_species(ptrPDMData, ptrPDMData->sv, species.index, species.coefficient, dlt_a_0);
				}

				// Update total propensity
				ptrPDMData->dTotalPropensity += dlt_a_0;
				sv.dTotalPropensity += dlt_a_0;
			}

			return true;
		}

		double UpdateModule_PDM::comb(int n, int coeff)
		{
			if (coeff == 0)	return 1.0;
			if (n <= 0)
			{
				n = 0;
			}
			int nt = n;
			for (int  i = 1 ; i < coeff; i++)
			{
				nt *= (n-i);
			}
			return nt;
		}
		
		void UpdateModule_PDM::update_structures_species(pssalib::datamodel::DataModel_PDM* ptrPDMData, UNSIGNED_INTEGER svix, UNSIGNED_INTEGER species, INTEGER species_coeff, REAL& dlt_a_0)
		{
			pssalib::datamodel::DataModel::Subvolume &sv = ptrPDMData->arSubvolumes[svix];
			pssalib::datamodel::DataModel_PDM::Subvolume_PDM &svpdm = ptrPDMData->arSubvolumes_PDM[svix];

			UNSIGNED_INTEGER U3_rowlen = ptrPDMData->arU3.get_cols(species);
			bool flag = true;

			// For each reaction that contains current species
			for(UNSIGNED_INTEGER l = 0; l < U3_rowlen; l++)
			{
				REAL change;
				const pssalib::datamodel::DataModel_PDM::PropensityIndex propIdx = ptrPDMData->arU3(species, l);
				UNSIGNED_INTEGER mu_af = ptrPDMData->aruL(propIdx.i, propIdx.j);
				
				int cN = 0;
				int reaction_coeff = 0;
				if (ptrPDMData->arUm.get_cols(mu_af) >= 2)
				{
					reaction_coeff = ptrPDMData->arUm(mu_af,1).coefficient;
					cN = sv.aruN[species];
				}
				else
				{
					reaction_coeff = (ptrPDMData->arUm(mu_af,0).coefficient < 0)? ptrPDMData->arUm(mu_af,0).coefficient + 1: ptrPDMData->arUm(mu_af,0).coefficient - 1;
					cN = sv.aruN[species]-1;
				}
				
				double comb1 = sv.ardC[mu_af]*comb(cN, abs(reaction_coeff))/factorial(abs(reaction_coeff));
				double comb2 = sv.ardC[mu_af]*comb(cN - species_coeff, abs(reaction_coeff))/factorial(abs(reaction_coeff));
				
/*				if(propIdx.i == species)
				{
					if (comb2/2 < svpdm.arPi(propIdx.i,propIdx.j)+0.5 && comb2/2 > svpdm.arPi(propIdx.i,propIdx.j))
					{
						std::cout << "NOOOOOOOOOOOO 3    " << comb2/2 << "   " << svpdm.arPi(propIdx.i,propIdx.j) << "\n\n\n";
						exit(-1);
					}
				}
				else
				{
					if (comb2 < svpdm.arPi(propIdx.i,propIdx.j)+0.5 && comb2 > svpdm.arPi(propIdx.i,propIdx.j))
					{
						std::cout << "NOOOOOOOOOOOOO 4 \n\n\n";
						exit(-1);
					}
				}*/
				
				// Calculate the change in reaction propensity
//				if (multi == true)
					change =  comb1 - comb2;
//				else
//					change = species_coeff * sv.ardC[mu_af];
				
				// If this is polymerization, account for it
				if(propIdx.i == species)
				{
					change /= 2.0;
					flag = false;
				}
				
				// Update the respective data structures
				svpdm.arPi(propIdx.i,propIdx.j) += change;
				
				if (svpdm.arPi(propIdx.i,propIdx.j) < 0.0 )
				{
					if (fabs(svpdm.arPi(propIdx.i,propIdx.j)) >= 1e-13)
					{
						std::cout << "Founded negative propensity value, simulation will stop\n";
						std::cout << svpdm.arPi << "   \n\n\n" << "change: " << change << "\n";
						exit(1);
					}
					else
					{
						svpdm.arPi(propIdx.i,propIdx.j) = 0.0;
					}
				}
				
				svpdm.arLambda[propIdx.i] += change;

				REAL temp = sv.aruN[propIdx.i] * svpdm.arLambda[propIdx.i];
				dlt_a_0 += temp - svpdm.arSigma[propIdx.i];
				svpdm.arSigma[propIdx.i] = temp;
			}

			// Update the group propensity of the affected species
			// (if it was not updated previously)
			if( flag )
			{
				REAL temp = sv.aruN[species] * svpdm.arLambda[species];
				dlt_a_0 += temp - svpdm.arSigma[species];
				svpdm.arSigma[species] = temp;
			}
		}
	}	// namespace update
}	// namespace pssalib
