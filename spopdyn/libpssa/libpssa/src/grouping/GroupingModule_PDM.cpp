/*
 * GroupingModule_PDM.cpp
 *
 *  Created on: 14 груд. 2010
 *      Author: sanyi
 */

#include "../../include/grouping/GroupingModule_PDM.h"
#include "../../include/update/UpdateModule_PDM.h"
#include "../../include/datamodel/DataModel_PDM.h"

namespace pssalib
{
	namespace grouping
	{
		GroupingModule_PDM::GroupingModule_PDM()
		{
		}

		GroupingModule_PDM::GroupingModule_PDM(GroupingModule &g) :
				GroupingModule(g)
		{
		}

		GroupingModule_PDM::~GroupingModule_PDM()
		{
		}

		////////////////////////////////
		// Grouping module methods

		//! Break-up non-elementary reactions?
		bool GroupingModule_PDM::do_breakup() const
		{
			return	false;
		}

		//! Initialize data structures (called before each trial)
		bool GroupingModule_PDM::initialize(pssalib::datamodel::DataModel* ptrData)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel_PDM* ptrPDMData = static_cast<pssalib::datamodel::DataModel_PDM*>(ptrData);

			// Call the baseclass method
			if(!GroupingModule::initialize(ptrData))
				return false;

			// Check for supported reactions
				
			for(UNSIGNED_INTEGER i = 0; i < ptrPDMData->unReactions; i++)
			{
				UNSIGNED_INTEGER nCols = ptrData->arUm.get_cols(i);
				
				if (nCols > 2)
				{
					std::cerr << "GroupingModule_PDM::initialize() - Error :" << " reaction " << i+1 << " has more than 2 reactants \n";

					bDataLoaded = false;
					return false;
				}
				
				if (nCols == 2 && abs(ptrData->arUm(i,0).coefficient) > 1)
				{
					std::cerr << "GroupingModule_PDM::initialize() - Error :" << " reaction " << i+1 << " has more than one spacies and needs a stoichiometry coefficent 1 for the first reactant \n";

					bDataLoaded = false;
					return false;
				}
			}
			
			// Clean the data members
			ptrPDMData->arU3.clear();
			ptrPDMData->aruL.clear();
			ptrPDMData->vQueuedReactions.clear();

			if(ptrPDMData->arSubvolumes_PDM != NULL)
				delete [] ptrPDMData->arSubvolumes_PDM;
			ptrPDMData->arSubvolumes_PDM = new pssalib::datamodel::DataModel_PDM::Subvolume_PDM[ptrData->unSubvolumes];
			
			// Preallocate memory
			UNSIGNED_INTEGER l = ptrPDMData->unReactions / ptrPDMData->unSpecies;
			if(0 == l) l = 1;
			ptrPDMData->arU3.reserve(ptrPDMData->unSpecies, l);
			ptrPDMData->aruL.reserve(ptrPDMData->unSpecies, l);
			
			///////////////////////////////////////////////////////
			// Initialize PDM data structures
			std::vector<pssalib::datamodel::DataModel::ChemicalSpecies>				idx_reactants;
			std::vector<pssalib::datamodel::DataModel::ChemicalSpecies>::iterator	minIt;
			pssalib::datamodel::DataModel_PDM::PropensityIndex 						idxPi;
			bool																	self_dependent;
			REAL 																	partprop;

			// Maximum number of reactants
			idx_reactants.reserve(3);

			// Partial propensity matrix indices
			for(UNSIGNED_INTEGER i = 0; i < ptrPDMData->unReactions; i++)
			{
				idx_reactants.clear();
				self_dependent = false;

				UNSIGNED_INTEGER nCols = ptrData->arUm.get_cols(i);
				
				// 
				
				for(UNSIGNED_INTEGER nu = 0; nu < nCols; nu++)
					idx_reactants.push_back(ptrData->arUm(i,nu));
				
				minIt = std::max_element(idx_reactants.begin(), idx_reactants.end());
				if(minIt->coefficient < -1)
				{
					self_dependent = true;
				}

				// Index of PI for the reaction
				idxPi.i = minIt->index;
				idxPi.j = ptrPDMData->aruL.get_cols(idxPi.i);

				// Account for self dependency
				if(self_dependent)
					ptrPDMData->arU3.push_back(minIt->index, idxPi);
				// Are any other species affected by the reaction?
				else if(idx_reactants.size() > 1)
				{
					idx_reactants.erase(minIt);
					// Account for them when changing propensity
					ptrPDMData->arU3.push_back(idx_reactants.begin()->index, idxPi);
				}
				ptrPDMData->aruL.push_back(idxPi.i, i);
			}

			// Diffusion reactions
			for(UNSIGNED_INTEGER nu = 1; nu < ptrData->unSpecies; nu++)
			{
				UNSIGNED_INTEGER reactionIdx = ptrPDMData->unReactions + nu - 1;
				ptrPDMData->aruL.push_back(nu, reactionIdx);
			}
			
			ptrPDMData->dTotalPropensity = 0.0;
			for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; ++si)
			{			  
				pssalib::datamodel::DataModel::Subvolume &sv = ptrPDMData->arSubvolumes[si];
				pssalib::datamodel::DataModel_PDM::Subvolume_PDM &svpdm = ptrPDMData->arSubvolumes_PDM[si];
				
				svpdm.arPi.reserve(ptrPDMData->unSpecies, l);

				// Partial propensity matrix
				for(UNSIGNED_INTEGER i = 0; i < ptrPDMData->unReactions; i++)
				{
					idx_reactants.clear();
					self_dependent = false;

					int tot = 0;
					UNSIGNED_INTEGER nCols = ptrData->arUm.get_cols(i);
					for(UNSIGNED_INTEGER nu = 0; nu < nCols; nu++)
					{idx_reactants.push_back(ptrData->arUm(i,nu));tot += abs(ptrData->arUm(i,nu).coefficient);}
					
					// Case 0: S1        ------>
					// Case 1: S1 +   S2 ------>
					// Case 2: S1 + X*S2 ------>
					
					if (idx_reactants.size() == 1)
						minIt = idx_reactants.begin();
					else if (idx_reactants.begin()->coefficient == (idx_reactants.begin() + 1)->coefficient)
						minIt = idx_reactants.begin()+1;
					else
						minIt = std::min_element(idx_reactants.begin(), idx_reactants.end());
					
					partprop = ptrPDMData->arSubvolumes[si].ardC[i];
					
					if (abs(minIt->coefficient) == tot)
					{
						// If reaction is of type      X*S1 ------> .......
					
						partprop *= ((REAL)/*(sv.aruN[maxIt->index]-1)*/pssalib::update::UpdateModule_PDM::comb(sv.aruN[minIt->index]-1,abs(minIt->coefficient)-1))/factorial(abs(minIt->coefficient));
						
						// Index of PI for the reaction
						
						idxPi.i = minIt->index;
					}
					else
					{
						// or                     S1 + X*S2 ------> .......
					
						partprop *= ((REAL)/*(sv.aruN[maxIt->index]-1)*/pssalib::update::UpdateModule_PDM::comb(sv.aruN[minIt->index],abs(minIt->coefficient)))/factorial(abs(minIt->coefficient));
						
						// Index of PI for the reaction
						
						idx_reactants.erase(minIt);
						idxPi.i = idx_reactants.begin()->index;
					}
					
					idxPi.j = svpdm.arPi.get_cols(idxPi.i);

					svpdm.arPi.push_back(idxPi.i, partprop);
				}
				
				// Diffusion reactions
				for(UNSIGNED_INTEGER nu = 0; nu < ptrData->unSpecies - 1; nu++)
				{
					UNSIGNED_INTEGER speciesIx = nu + 1;
					
					REAL lumpedP = sv.ardD[nu];
					lumpedP *= sv.unDiffusionReactions;
					svpdm.arPi.push_back(speciesIx, lumpedP);
				}

				// Calculate Sigma and Lambda
				svpdm.arSigma = new REAL[ptrPDMData->unSpecies];
				svpdm.arLambda = new REAL[ptrPDMData->unSpecies];

				memset(svpdm.arSigma, 0, sizeof(REAL) * ptrPDMData->unSpecies);
				memset(svpdm.arLambda, 0, sizeof(REAL) * ptrPDMData->unSpecies);

				sv.dTotalPropensity = 0.0;
				for(UNSIGNED_INTEGER i = 0; i < ptrPDMData->unSpecies; i++)
				{
					UNSIGNED_INTEGER n = svpdm.arPi.get_cols(i);
					for(UNSIGNED_INTEGER j = 0; j < n; j++)
						svpdm.arLambda[i] += svpdm.arPi(i,j);
					svpdm.arSigma[i] = svpdm.arLambda[i] * ptrPDMData->arSubvolumes[si].aruN[i];
					sv.dTotalPropensity += ptrPDMData->arSubvolumes_PDM[si].arSigma[i];
					ptrPDMData->dTotalPropensity += ptrPDMData->arSubvolumes_PDM[si].arSigma[i];
				}
			}

			return true;
		}
	}
}
