/*
 * GroupingModule.cpp
 *
 *  Created on: 14 груд. 2010
 *      Author: sanyi
 */

#include "../../include/grouping/GroupingModule.h"
#include "../../include/datamodel/SimulationInfo.h"

#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>



#define	BREAKUP_REACTION_K		((REAL)100000.0)
#define	BREAKUP_REACTION_INC	((REAL)10.0)


LIBSBML_CPP_NAMESPACE_USE

namespace pssalib
{
	namespace grouping
	{
	  
		// This is an iterator help class, suppose you have this reaction
		// 2X(1) + 3X(2) + 1X(3)
		// while (hasNext()) {getNextReactant()}
		// return the sequence X(1) + X(1) + X(2) + X(2) + X(2) + X(3)
		// 0 stochiometry coefficents are not allowed
	  
		class ReactionIt
		{
			private:
		  
				Reaction * r;
				int actualSp;
				int actualSto;
				SpeciesReference * ref;
				bool bReverseLoop;
				int nsp;
				
			public:
		  
				ReactionIt(Reaction * rt, bool bRL)
				{
					r = rt;
					actualSp = 0;
					bReverseLoop = bRL;
					
					ref = bReverseLoop ? r->getProduct(0) : r->getReactant(0);
					nsp = bReverseLoop ? r->getNumProducts() : r->getNumReactants();
					actualSto = ref->getStoichiometry();
				};
				
				SpeciesReference * getNextReactant()
				{
					if (actualSto > 0)
					{
						actualSto--;
						return ref;
					}
					else
					{
						actualSp++;
						if (actualSp < nsp)
						{
							ref = r->getReactant(actualSp);
							actualSto = ref->getStoichiometry();
							actualSto--;
							return ref;
						}
					}
				};
				
				bool hasNext()
				{
					if (actualSto > 0)
						return true;
					else
					{
						return (actualSp+1 < nsp);
					}
				}
			
		};
	  
	  
		// Constructor
		GroupingModule::GroupingModule()
		{
		}

		// Copy constructor
		GroupingModule::GroupingModule(GroupingModule &g) :
				bDataLoaded(g.bDataLoaded),
				mapSp2Idx(g.mapSp2Idx),
				arIdx2Sp(g.arIdx2Sp)
		{
		}

		// Destructor
		GroupingModule::~GroupingModule()
		{
		}
		
		///////////////////////////////////////////////////////////////////////
		// Grouping module methods
		
		void GroupingModule::printsp(int coef, std::string name)
		{
			if (coef != 1)
				std::cout << abs(coef) << name;
			else
				std::cout << name;
		}
		
		void GroupingModule::PrintSimulationInformation(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			std::cout << " ----------- Reactions ------------- \n\n";
			
			pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();
			
			for (int i = 0 ; i < ptrData->arUp.get_rows() ; i++)
			{
				std::cout << i+1 << ")    ";
			  
				if (ptrData->arUm.get_cols(i) != 0)
					printsp(abs(ptrData->arUm(i,0).coefficient),ptrData->speciesName[ptrData->arUm(i,0).index]);

				for (int j = 1 ; j  < ptrData->arUm.get_cols(i) ; j++)
				{
					std::cout << " + ";
					printsp(abs(ptrData->arUm(i,j).coefficient),ptrData->speciesName[ptrData->arUm(i,j).index]);
				}
				
				std::cout << "   ------>(" << ptrData->arSubvolumes->ardC[i] << ")   ";
				
				if (ptrData->arUp.get_cols(i) != 0)
					printsp(abs(ptrData->arUp(i,0).coefficient),ptrData->speciesName[ptrData->arUp(i,0).index]);
				
				for (int j = 1 ; j  < ptrData->arUp.get_cols(i) ; j++)
				{
					std::cout << " + ";
					printsp(abs(ptrData->arUp(i,j).coefficient),ptrData->speciesName[ptrData->arUp(i,j).index]);
				}
				
				std::cout << "\n";
			}
			
			std::cout << " ----------- Diffusion ------------- \n\n";
			
			for (int i = 1 ; i < ptrData->unSpecies ; i++)
			{
				std::cout << i << ")    ";
			  
				printsp(1,ptrData->speciesName[i]);
				
				std::cout << "   ------>(" << ptrData->arSubvolumes->ardD[i-1] << ")   ";
				
				printsp(1,ptrData->speciesName[i]);
				
				std::cout << "\n";
			}
		}
		
		bool GroupingModule::preinitialize(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();

			if(NULL == ptrData||!ptrSimInfo->isDataLoaded())
				throw std::runtime_error("GroupingModule::read_SBML() - data structures not initialized!");

			// Clean-up previous model
			ptrData->Cleanup();
			ptrData->arU.free();
			ptrData->arUp.free();
			ptrData->arUm.free();

			// Get the model
			Model *model = ptrSimInfo->getSBMLModel();

			REAL volume = (REAL)model->getCompartment(0)->getVolume();

			std::vector<ReactionInfo> vReactionInfo;
			if (!parse_parameters(ptrSimInfo, model, vReactionInfo, volume))
				return false;

			if (ptrSimInfo->bVerbose)
				std::cout << "GroupingModule::preinitialize() : Done with parameters" << std::endl;

			ptrData->unSpecies	 = model->getNumSpecies() + 1 + unNumAdditionalTotal;
			ptrData->unReactions = model->getNumReactions() + 2 * unNumAdditionalTotal + unReverseReactions;

			ptrData->arU.reserve(ptrData->unReactions, 2);
			ptrData->arUp.reserve(ptrData->unReactions, 1);
			ptrData->arUm.reserve(ptrData->unReactions, 1);
			
			ptrData->arU_r.reserve(ptrData->unReactions, 2);
			ptrData->arUp_r.reserve(ptrData->unReactions, 1);
			ptrData->arUm_r.reserve(ptrData->unReactions, 1);

			switch(ptrSimInfo->eVolumeType)
			{
			case datamodel::SimulationInfo::VT_SingleSubvolume:
				if (!setup_volumes_single(ptrData, ptrSimInfo->eBoundaryConditions, ptrSimInfo->dOmega))
					return false;
				break;

			case datamodel::SimulationInfo::VT_Homogeneous2d:
				if (!setup_volumes_homogeneous_2d(ptrData, ptrSimInfo->eBoundaryConditions, ptrSimInfo->unNumGridPoints, ptrSimInfo->dOmega))
					return false;
				break;

			case datamodel::SimulationInfo::VT_Homogeneous3d:
				if (!setup_volumes_homogeneous_3d(ptrData, ptrSimInfo->eBoundaryConditions, ptrSimInfo->unNumGridPoints, ptrSimInfo->dOmega))
					return false;
				break;
				  
			default:
				std::cerr << "GroupingModule::preinitialize() : Invalid volume type: " << ptrSimInfo->eVolumeType << std::endl;
				return false;
			}

			if (!initialize_mappings(ptrSimInfo, ptrData, model))
				return false;
			if (!fill_mappings(ptrSimInfo, ptrData, model, vReactionInfo))
				return false;
			if (!fill_diffusion(ptrSimInfo))
				return false;
			if (!adjust_reaction_rates_by_volume(ptrData))
				return false;

			bDataLoaded = true;

			if (ptrSimInfo->bVerbose)
			{
				std::cout << "arU  :" << std::endl << ptrData->arU  << std::endl;
				std::cout << "arUm :" << std::endl << ptrData->arUm << std::endl;
				std::cout << "arUp :" << std::endl << ptrData->arUp << std::endl;
			}

			PrintSimulationInformation(ptrSimInfo);
			
			return true;
		}

		//! Initialize data structures (called before each trial)
		bool GroupingModule::initialize(pssalib::datamodel::DataModel* ptrData)
		{
			for (UNSIGNED_INTEGER i = 0; i < ptrData->unSubvolumes; i++)
			{
				ptrData->arSubvolumes[i].aruN[0] = 1;
				ptrData->arSubvolumes[i].dTotalPropensity = 0.0;
			}
			ptrData->dTotalPropensity = 0.0;

			return true;
		}

		bool GroupingModule::postinitialize(pssalib::datamodel::DataModel* ptrData)
		{
			// Calculate minimum value
			ptrData->crsdVolume.bins.resize(ptrData->unSubvolumes);
			ptrData->crsdVolume.minValue = std::numeric_limits<REAL>::max();
			for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; ++si)
			{
				pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[si];

				REAL minPropensity = std::numeric_limits<REAL>::max();
				for(UNSIGNED_INTEGER i = 0; i < ptrData->unReactions; i++)
				{
					REAL temp = sv.ardC[i];
					
					// zero reaction continue
					
					if (sv.ardC[i] == 0)
						continue;
					
					if (ptrData->arUp.get_cols(i) == 2) {
						if (ptrData->arUp(i, 0) == ptrData->arUp(i, 1)) {
							temp *= 0.5;
						}
					}

					if (temp < minPropensity)
						minPropensity = temp;
				}
				
				for(UNSIGNED_INTEGER nu = 0; nu < ptrData->unSpecies-1; nu++)
				{
					REAL temp = sv.ardD[nu] / (sv.dH * sv.dH);
					if (temp < minPropensity)
						minPropensity = temp;
				}

				// Don't need to add diffusion; all diffusion reactions have a minimum propensity of 0.

				if (minPropensity < ptrData->crsdVolume.minValue)
					ptrData->crsdVolume.minValue = minPropensity;
			}
			
			// Distribute the subvolume propensities into the bins.
			for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; ++si)
			{
				UNSIGNED_INTEGER k;
				pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[si];
				k = (UNSIGNED_INTEGER)std::floor(fabs(log2(sv.dTotalPropensity / ptrData->crsdVolume.minValue))) + 1;
				ptrData->crsdVolume.updateValue(k, si, sv.dTotalPropensity);
			}

			return true;
		}
		
		double GroupingModule::eval(const ASTNode * nd, double v1, double v2)
		{
			switch (nd->getType())
			{
				case AST_PLUS:
					return v1 + v2;
					break;
				case AST_MINUS:
					return v1 - v2;
					break;
				case AST_TIMES:
					return v1 * v2;
					break;
				case AST_DIVIDE:
					return v1 / v2;
					break;
				case AST_FUNCTION_POWER:
				case AST_POWER:
					return pow(v1,v2);
					break;
				default:
					std::cerr << "GroupingModule::eval - Error : Unrecognized AST_NODE" << std::endl;
					eval_error = true;
			}
		}
		
		// Search inside the parameter search
		
		double GroupingModule::searchPar(std::string str, KineticLaw * L, Model * model)
		{
			for (int i = 0 ; i < model->getNumParameters() ; i++)
			{
				Parameter * par = model->getParameter(i);
				
				std::string name = par->getName();
				
				if (par->getName().compare(str) == 0)
					return par->getValue();
			}
			
			for (int i = 0 ; i < L->getNumParameters() ; i++)
			{
				Parameter * par = L->getParameter(i);
				
				if (par->getName().compare(str) == 0)
					return par->getValue();
			}

			std::cerr	<< "GroupingModule::searchPar - Error : Unknown parameter " << str << std::endl;
			eval_error = true;
			
			return 0.0;
		}
		
		double GroupingModule::evaluateMath(const ASTNode * Mnode, KineticLaw * L, Model * model)
		{
			int nop = 0;
			double val = 0.0;
			double tval = 0.0;
			for (unsigned int i = 0 ; i < Mnode->getNumChildren() ; i++)
			{
				ASTNode * Mtnode = Mnode->getChild(i);
				int test = Mtnode->getType();
				if (Mtnode->getType() == AST_PLUS   || Mtnode->getType() == AST_MINUS ||
				    Mtnode->getType() == AST_MINUS  || Mtnode->getType() == AST_TIMES ||
				    Mtnode->getType() == AST_DIVIDE || Mtnode->getType() == AST_POWER ||
				    Mtnode->getType() == AST_FUNCTION_POWER )
				{
					if (nop == 0)
					{
						val = evaluateMath(Mtnode,L,model);
						nop++;
					}
					else
					{
						tval = evaluateMath(Mtnode,L,model);
						val = eval(Mnode,tval,val);
					}
				}
				else if (Mtnode->getType() == AST_NAME)
				{
					if (nop == 0)
					{
						val = searchPar(Mtnode->getName(),L,model);
						nop++;
					}
					else
					{
						tval = searchPar(Mtnode->getName(),L,model);
						val = eval(Mnode,val,tval);
					}
				}
				else if (Mtnode->getType() == AST_INTEGER)
				{
					if (nop == 0)
					{
						val = (double)Mtnode->getInteger();
						nop++;
					}
					else
					{
						tval = (double)Mtnode->getInteger();
						val = eval(Mnode,val,tval);
					}
				}
				else if (Mtnode->getType() == AST_REAL)
				{
					if (nop == 0)
					{
						val = (double)Mtnode->getReal();
						nop++;
					}
					else
					{
						tval = (double)Mtnode->getReal();
						val = eval(Mnode,val,tval);
					}
				}
				else
				{
					std::cerr	<< "GroupingModule::evaluateMath - Error : Unrecognized AST_NODE    " << Mtnode->getType() << std::endl;
					eval_error = true;
				}
			}
			
			return val;
		}
		
		bool GroupingModule::parse_parameters(pssalib::datamodel::SimulationInfo* ptrSimInfo, Model* model, std::vector<ReactionInfo>& vReactionInfo, REAL volume)
		{
			Reaction 								*reaction;
			SpeciesReference 						*speciesReference;
			KineticLaw 							*kineticLaw;
			Parameter								*parameter;

			
			unNumAdditionalTotal = 0;
			unReverseReactions = 0;

			bool bReverseLoop = false;

			{
			  ReactionInfo		ri;
			  
			  
				bool 				rateSet = false, cmuSet = false;
				REAL				dReactRate = 0.0, dRevReactRate = 0.0;
				for (UNSIGNED_INTEGER mu = 0, n = 0, unNumAdditional = 0, unParams = 0, unL = 0; mu < model->getNumReactions(); ++mu, ++unL)
				{
					unNumAdditional = 0;				reaction = model->getReaction(mu);
					n = bReverseLoop ? reaction->getNumProducts() : reaction->getNumReactants();

					/// GET Reaction annotation
					if (reaction->isSetAnnotation()){
					  XMLNode * annotation = reaction->getAnnotation();
					  ri.spatial_c = getSpatialC(annotation,ptrSimInfo->unNumGridPoints);
						}
					else{
					  ri.spatial_c = NULL;
						}
					//// 

					
					if (ptrSimInfo->bVerbose)
					{
						std::cout	<< "VERBOSE : reaction " << reaction->getId() << (bReverseLoop ? "(reverse)" : "") << " is"
									<< (reaction->getReversible() ? " " : " not ") << "reversible , loop " << unL << std::endl;
					}

					if(!bReverseLoop)
					{
						rateSet = false;	cmuSet = false;	ri.clear();
						kineticLaw = reaction->getKineticLaw();	unParams = kineticLaw->getNumParameters();

						dReactRate = 0.0, dRevReactRate = 0.0;
						 
						if (kineticLaw->isSetMath() == true)
						{
							// Get Math
						  
							const ASTNode * Mnodes = kineticLaw->getMath();

							eval_error = false;
							dReactRate = evaluateMath(Mnodes,kineticLaw,model);
							
							// Math evaluation error
							
							if (eval_error == true)		return false;
							
							rateSet = true;
						}
						
						for (UNSIGNED_INTEGER i = 0; i < unParams; ++i)
						{
							parameter = kineticLaw->getParameter(i);

							if(parameter->getName().compare("c") == 0)	// rate constant
							{
								ri.c	= ((REAL)parameter->getValue());
					
								if (ptrSimInfo->bVerbose)
									std::cout << "\tc = " << ri.c << std::endl;

								cmuSet = true;
							}
							else if(parameter->getName().compare("k") == 0)	// rate constant
							{
								dReactRate = ((REAL)parameter->getValue());

								if (ptrSimInfo->bVerbose)
									std::cout << "\tk = " << dReactRate << std::endl;

								rateSet = true;
							}
							else if(parameter->getName().compare("tau") == 0)	// time delay
							{
								// FIXME Implement a way to distinguish between D1 & D2 classes
								ri.rc	= pssalib::datamodel::DataModel::D2;
								// Set delay
								ri.d	= ((REAL)parameter->getValue());

								if (ptrSimInfo->bVerbose)
									std::cout << "\tdelay = " << ri.d << std::endl;
							}
							// Only if the reaction is reversible
							if(reaction->getReversible())
							{
								if(parameter->getName().compare("c_r") == 0)	// rate constant
								{
									ri.c_r = ((REAL)parameter->getValue());
									
									if (ptrSimInfo->bVerbose)
										std::cout << "\tc(reverse) = " << ri.d << std::endl;
								}
								else if(parameter->getName().compare("k_r") == 0)	// rate constant
								{
									dRevReactRate = ((REAL)parameter->getValue());
									if (ptrSimInfo->bVerbose)
										std::cout << "\tk(reverse) = " << dRevReactRate << std::endl;
								}
							}
						}

						// Check if any reaction parameters are set
						if((!rateSet)&&(!cmuSet))
						{
							std::stringstream	ssTemp;
							ssTemp	<< "GroupingModule::preinitialize() - Error :" << " reaction " << reaction->getId()
									<< " should have a rate constant associated with it.";
							ptrSimInfo->reportError(ssTemp.str().c_str());

							bDataLoaded = false;
							return false;
						}
						else if(rateSet&&cmuSet)
						{
							std::stringstream	ssTemp;
							ssTemp	<< "GroupingModule::preinitialize() - Error :" << "reaction " << reaction->getId()
									<< " has both a rate constant and a specific probability "
									<< "rate associated with it." << std::endl;
							ptrSimInfo->reportError(ssTemp.str().c_str());
							bDataLoaded = false;
							return false;
						}

						// Copy the reaction rate for the reversible reaction if not set
						if(0.0 == dRevReactRate)
							dRevReactRate	= dReactRate;
						// Copy the specific probability rate for the reversible reaction if not set
						if(0.0 == ri.c_r)
							ri.c_r			= ri.c;
					}

					if(rateSet)
					{
/*						UNSIGNED_INTEGER	unNuProd = 1;
						for(UNSIGNED_INTEGER j = 0, nu = 0; j < n; j++)
						{
							speciesReference	 = bReverseLoop ? reaction->getProduct(j) : reaction->getReactant(j);
							nu					 = 1 * std::abs(speciesReference->getStoichiometry());
							unNuProd 			*= factorial(nu);
							unNumAdditional		+= nu;
						}

						if (ptrSimInfo->bVerbose)
							std::cout << "\tconverting reaction rate into specific probability rate..." << std::endl;
*/
						if(bReverseLoop)
							ri.c_r	= /*( */dRevReactRate /* * (REAL)unNuProd) / pow(volume, unNumAdditional - 1.0)*/;
						else
							ri.c	= /*(*/dReactRate  /* * (REAL)unNuProd) / pow(volume, unNumAdditional - 1.0)*/;
					}
					else
					{
						for(UNSIGNED_INTEGER j = 0, nu = 0; j < n; j++)
						{
							speciesReference	 = bReverseLoop ? reaction->getProduct(j) : reaction->getReactant(j);
							nu					 = 1 * std::abs(speciesReference->getStoichiometry());
							unNumAdditional		+= nu;
						}
					}

					if(do_breakup()&&(unNumAdditional > 2))
					{
						if(bReverseLoop)
							ri.n_r	= unNumAdditional - 2;
						else
							ri.n	= unNumAdditional - 2;
						unNumAdditionalTotal	+= unNumAdditional - 2;
					}

					if(!bReverseLoop&&reaction->getReversible())
					{
						bReverseLoop = true;	++unReverseReactions;	--mu;	continue;
					}
					else if(bReverseLoop)
						bReverseLoop = false;

					vReactionInfo.push_back(ri);
				}
			}

			return true;
		}
		
		bool GroupingModule::setup_volumes_single(pssalib::datamodel::DataModel* ptrData, pssalib::datamodel::SimulationInfo::BoundaryConditionsType bc, REAL omega)
		{
			ptrData->unSubvolumes = 1;
			ptrData->unMiddleVolumeIx = 0;
			REAL svomega = omega;
			REAL h = omega;
			return initialize_volumes(ptrData, 0, svomega, h);
		}
		
		bool GroupingModule::setup_volumes_homogeneous_2d(pssalib::datamodel::DataModel* ptrData, pssalib::datamodel::SimulationInfo::BoundaryConditionsType bc, UNSIGNED_INTEGER numGridPoints, REAL omega)
		{
			ptrData->unSubvolumes = numGridPoints * numGridPoints;
			ptrData->unMiddleVolumeIx = (numGridPoints / 2) * numGridPoints + numGridPoints / 2;
			REAL svomega = omega / ptrData->unSubvolumes;
			REAL h = sqrt(svomega);

			if (!initialize_volumes(ptrData, 4, svomega, h))
				return false;
			
			if (bc == datamodel::SimulationInfo::BC_Periodic)
			{
				// Periodic boundary conditions.
				for (UNSIGNED_INTEGER y = 0; y < numGridPoints; y++)
				{
					for (UNSIGNED_INTEGER x = 0; x < numGridPoints; x++)
					{
						UNSIGNED_INTEGER index = y * numGridPoints + x;
						pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[index];
					
						UNSIGNED_INTEGER x1 = ((x-1) + numGridPoints) % numGridPoints;
						UNSIGNED_INTEGER x2 = ((x+1) + numGridPoints) % numGridPoints;
						UNSIGNED_INTEGER y1 = ((y-1) + numGridPoints) % numGridPoints;
						UNSIGNED_INTEGER y2 = ((y+1) + numGridPoints) % numGridPoints;

						sv.arNeighborVolumes[0] = y1 * numGridPoints + x;
						sv.arNeighborVolumes[1] = y2 * numGridPoints + x;
						sv.arNeighborVolumes[2] = y * numGridPoints + x1;
						sv.arNeighborVolumes[3] = y * numGridPoints + x2;
					}
				}
			}
			else if (bc == datamodel::SimulationInfo::BC_Reflexive)
			{
				// Reflexive boundary conditions.
				for (UNSIGNED_INTEGER y = 0; y < numGridPoints; y++)
				{
					for (UNSIGNED_INTEGER x = 0; x < numGridPoints; x++)
					{
						UNSIGNED_INTEGER index = y * numGridPoints + x;
						pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[index];
					
						UNSIGNED_INTEGER x1 = (std::max)(0, (int)x-1);
						UNSIGNED_INTEGER x2 = (std::min)(x+1, numGridPoints-1);
						UNSIGNED_INTEGER y1 = (std::max)(0, (int)y-1);
						UNSIGNED_INTEGER y2 = (std::min)(y+1, numGridPoints-1);

						sv.arNeighborVolumes[0] = y1 * numGridPoints + x;
						sv.arNeighborVolumes[1] = y2 * numGridPoints + x;
						sv.arNeighborVolumes[2] = y * numGridPoints + x1;
						sv.arNeighborVolumes[3] = y * numGridPoints + x2;
					}
				}
			}
			else
			{
				std::cerr << "GroupingModule::setup_volumes_homogeneous_2d: Invalid boundary conditions: " << bc << std::endl;
				return false;
			}

			return true;
		}
		
		bool GroupingModule::setup_volumes_homogeneous_3d(pssalib::datamodel::DataModel* ptrData, pssalib::datamodel::SimulationInfo::BoundaryConditionsType bc, UNSIGNED_INTEGER numGridPoints, REAL omega)
		{
			ptrData->unSubvolumes = numGridPoints * numGridPoints * numGridPoints;
			ptrData->unMiddleVolumeIx = (numGridPoints / 2) * numGridPoints * numGridPoints + (numGridPoints / 2) * numGridPoints + numGridPoints / 2;
			REAL svomega = omega / ptrData->unSubvolumes;
			REAL h = pow(svomega, 1./3.);

			if (!initialize_volumes(ptrData, 6, svomega, h))
				return false;
			
			if (bc == datamodel::SimulationInfo::BC_Periodic)
			{
				// Periodic boundary conditions.
				for (UNSIGNED_INTEGER z = 0; z < numGridPoints; z++)
				{
					for (UNSIGNED_INTEGER y = 0; y < numGridPoints; y++)
					{
						for (UNSIGNED_INTEGER x = 0; x < numGridPoints; x++)
						{
							UNSIGNED_INTEGER index = z * numGridPoints * numGridPoints + y * numGridPoints + x;
							pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[index];
					
							UNSIGNED_INTEGER x1 = ((x-1) + numGridPoints) % numGridPoints;
							UNSIGNED_INTEGER x2 = ((x+1) + numGridPoints) % numGridPoints;
							UNSIGNED_INTEGER y1 = ((y-1) + numGridPoints) % numGridPoints;
							UNSIGNED_INTEGER y2 = ((y+1) + numGridPoints) % numGridPoints;
							UNSIGNED_INTEGER z1 = ((z-1) + numGridPoints) % numGridPoints;
							UNSIGNED_INTEGER z2 = ((z+1) + numGridPoints) % numGridPoints;

							sv.arNeighborVolumes[0] = z1 * numGridPoints * numGridPoints + y * numGridPoints + x;
							sv.arNeighborVolumes[1] = z2 * numGridPoints * numGridPoints + y * numGridPoints + x;
							sv.arNeighborVolumes[2] = z * numGridPoints * numGridPoints + y1 * numGridPoints + x;
							sv.arNeighborVolumes[3] = z * numGridPoints * numGridPoints + y2 * numGridPoints + x;
							sv.arNeighborVolumes[4] = z * numGridPoints * numGridPoints + y * numGridPoints + x1;
							sv.arNeighborVolumes[5] = z * numGridPoints * numGridPoints + y * numGridPoints + x2;
						}
					}
				}
			}
			else if (bc == datamodel::SimulationInfo::BC_Reflexive)
			{
				// Reflexive boundary conditions.
				for (UNSIGNED_INTEGER z = 0; z < numGridPoints; z++)
				{
					for (UNSIGNED_INTEGER y = 0; y < numGridPoints; y++)
					{
						for (UNSIGNED_INTEGER x = 0; x < numGridPoints; x++)
						{
							UNSIGNED_INTEGER index = z * numGridPoints * numGridPoints + y * numGridPoints + x;
							pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[index];
					
							UNSIGNED_INTEGER x1 = (std::max)(0, (int)x-1);
							UNSIGNED_INTEGER x2 = (std::min)(x+1, numGridPoints-1);
							UNSIGNED_INTEGER y1 = (std::max)(0, (int)y-1);
							UNSIGNED_INTEGER y2 = (std::min)(y+1, numGridPoints-1);
							UNSIGNED_INTEGER z1 = (std::max)(0, (int)z-1);
							UNSIGNED_INTEGER z2 = (std::min)(z+1, numGridPoints-1);

							sv.arNeighborVolumes[0] = z1 * numGridPoints * numGridPoints + y * numGridPoints + x;
							sv.arNeighborVolumes[1] = z2 * numGridPoints * numGridPoints + y * numGridPoints + x;
							sv.arNeighborVolumes[2] = z * numGridPoints * numGridPoints + y1 * numGridPoints + x;
							sv.arNeighborVolumes[3] = z * numGridPoints * numGridPoints + y2 * numGridPoints + x;
							sv.arNeighborVolumes[4] = z * numGridPoints * numGridPoints + y * numGridPoints + x1;
							sv.arNeighborVolumes[5] = z * numGridPoints * numGridPoints + y * numGridPoints + x2;
						}
					}
				}
			}
			else
			{
				std::cerr << "GroupingModule::setup_volumes_homogeneous_2d: Invalid boundary conditions: " << bc << std::endl;
				return false;
			}

			return true;
		}

		bool GroupingModule::initialize_volumes(pssalib::datamodel::DataModel* ptrData, UNSIGNED_INTEGER numDiffusionReactions, REAL omega, REAL h)
		{
			ptrData->arSubvolumes = new pssalib::datamodel::DataModel::Subvolume[ptrData->unSubvolumes];

			for (UNSIGNED_INTEGER i = 0; i < ptrData->unSubvolumes; i++)
			{
				ptrData->arSubvolumes[i].dOmega					= omega;
				ptrData->arSubvolumes[i].dH						= h;
				ptrData->arSubvolumes[i].unDiffusionReactions	= numDiffusionReactions;
				ptrData->arSubvolumes[i].aruN 					= new INTEGER[ptrData->unSpecies];
				ptrData->arSubvolumes[i].ardC					= new REAL[ptrData->unReactions];
				ptrData->arSubvolumes[i].ardD					= new REAL[ptrData->unSpecies-1];
				ptrData->arSubvolumes[i].arReactClass			= new pssalib::datamodel::DataModel::ReactionClass[ptrData->unReactions + ptrData->unSpecies - 1];
				
				memset(ptrData->arSubvolumes[i].aruN, 0, sizeof(INTEGER) * ptrData->unSpecies);
				memset(ptrData->arSubvolumes[i].ardC, 0, sizeof(REAL) * ptrData->unReactions);
				memset(ptrData->arSubvolumes[i].ardD, 0, sizeof(REAL) * (ptrData->unSpecies-1));
				memset(ptrData->arSubvolumes[i].arReactClass, 0, sizeof(pssalib::datamodel::DataModel::ReactionClass) * (ptrData->unReactions + ptrData->unSpecies - 1));
				
				if (numDiffusionReactions > 0)
				{
					ptrData->arSubvolumes[i].arNeighborVolumes = new UNSIGNED_INTEGER[numDiffusionReactions];
					memset(ptrData->arSubvolumes[i].arNeighborVolumes, 0, sizeof(UNSIGNED_INTEGER) * numDiffusionReactions);
				}
			}

			return true;
		}

		bool GroupingModule::initialize_mappings(pssalib::datamodel::SimulationInfo* ptrSimInfo, pssalib::datamodel::DataModel* ptrData, Model* model)
		{
			arIdx2Sp.clear();
			mapSp2Idx.clear();
			arIdx2Sp.reserve(model->getNumSpecies());
			mapSp2Idx.rehash(model->getNumSpecies());
			ptrData->speciesName.push_back(std::string("0"));
			for(UNSIGNED_INTEGER i = 0; i < model->getNumSpecies(); ++i)
			{
			  Species *species = model->getSpecies(i);
				
				if (ptrSimInfo->bVerbose)
					std::cout << "Species : id = " << species->getId() << " : name = '" << species->getName() << "' : init_amount = " << species->getInitialAmount() << std::endl;

				mapSp2Idx.insert( MAP_SP2IDX::value_type(species->getId(), i + 1) );
				arIdx2Sp.push_back(species->getId());
				std::string name(species->getName());
				ptrData->speciesName.push_back(name);
				
				if (ptrSimInfo->eInitialPop == pssalib::datamodel::SimulationInfo::IP_Distribute) {
					for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; si++) {
						ptrData->arSubvolumes[si].aruN[i+1] = ((INTEGER)species->getInitialAmount()) / ptrData->unSubvolumes;
					}
				} else if (ptrSimInfo->eInitialPop == pssalib::datamodel::SimulationInfo::IP_Concentrate) {
					if (ptrData->unMiddleVolumeIx >= 0 && ptrData->unMiddleVolumeIx < ptrData->unSubvolumes) {
						ptrData->arSubvolumes[ptrData->unMiddleVolumeIx].aruN[i+1] = (INTEGER)species->getInitialAmount();
					} else {
						std::cerr << "GroupingModule::initialize_mappings: Invalid middle volume index" << std::endl;
					}
				} else if (ptrSimInfo->eInitialPop == pssalib::datamodel::SimulationInfo::IP_Multiply) {
					for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; si++) {
						ptrData->arSubvolumes[si].aruN[i+1] = (INTEGER)species->getInitialAmount();
					}
				} else if (ptrSimInfo->eInitialPop == pssalib::datamodel::SimulationInfo::IP_UserDefined) {
				  
					if (ptrSimInfo->eInitialUserDefined == NULL) {
						std::cerr << "GroupingModule::initialize_mappings: No initial population callback" << std::endl;
					}
					pssalib::datamodel::DataModel::SubvolumesInit * in = new pssalib::datamodel::DataModel::SubvolumesInit(ptrData->unSpecies,ptrData->unSubvolumes,ptrData->arSubvolumes);
					
					try
					{
						ptrSimInfo->eInitialUserDefined(in);
					}
					catch (std::string exception)
					{
						std::cerr << exception;
					}
					
					delete(in);
				}
			}

			return true;
		}

		bool GroupingModule::fill_mappings(pssalib::datamodel::SimulationInfo* ptrSimInfo, pssalib::datamodel::DataModel* ptrData, Model* model, const std::vector<ReactionInfo>& vReactionInfo)
		{
			bool bReverseLoop = false;

			pssalib::datamodel::DataModel::ChemicalSpecies cs;
			for (UNSIGNED_INTEGER mu = 0, mu_actual = 0, unEntries = 0, unIdxAdditionalSpecies = model->getNumSpecies() + 1; mu < model->getNumReactions(); ++mu, ++mu_actual)
			{
				Reaction* reaction	= model->getReaction(mu);
				const ReactionInfo& ri = vReactionInfo[mu];

				unEntries = bReverseLoop ? reaction->getNumProducts() : reaction->getNumReactants();
				if(unEntries == 0)
				{
					std::cerr << "Warning : reaction " << reaction->getId() << (bReverseLoop ? "(reverse)" : "") << " has no reactants! Assuming 'reservoir'." << std::endl;
					cs.index = 0;
					cs.coefficient = -1;
					ptrData->arU.push_back(mu_actual, cs);
					ptrData->arUm.push_back(mu_actual, cs);

					for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; si++)
					{
						pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[si];
						sv.ardC[mu_actual]			= (bReverseLoop ? ri.c_r : ri.c);
						sv.arReactClass[mu_actual]	= pssalib::datamodel::DataModel::D0;
					}
				}
				else
				{
					if (ptrSimInfo->bVerbose)
						std::cout << "Reaction " << reaction->getId() << " has " << reaction->getNumReactants() << " reactants : " << std::endl << '\t';

					// Check if a break-up is needed
					if(do_breakup()&&((!bReverseLoop&&(ri.n > 0))||(bReverseLoop&&(ri.n_r > 0)))) // do break-up
					{
						// Get the first two species of the first bimolecular reaction

						UNSIGNED_INTEGER	unNumAdditional = 0, unSpIdx = 0;
						
						//
						// The breakdown of a n reaction into bimolecular reaction
						// require that you have a fast system and a slow system
						// such that you can solve the Z species (the fast system)
						// concentration based on the state of the system. This mean
						// that every time you fire some other reaction (you change)
						// the state of the system your Z species it should quickly termalize
						// accordly to the new state of the system. This mean that your
						// have to fire a lot of reaction of your fast system, reducing
						// drammaticaly the efficency of your simulation
						//
						
						//
						//   Search a suitable max rate
						//   double max rate
						//
						
						REAL				C	= pow( ((REAL)BREAKUP_REACTION_K) * (bReverseLoop ? ri.c_r : ri.c),
								((REAL)(bReverseLoop ? ri.n_r : ri.n)) / ((REAL)((bReverseLoop ? ri.n_r : ri.n) + 1)) );

						ReactionIt reactionIt(reaction,bReverseLoop);
						
						SpeciesReference *speciesReference1 = reactionIt.getNextReactant();
						
						// Find the index of the first reactant
						MAP_SP2IDX::iterator	it1 = mapSp2Idx.find(speciesReference1->getSpecies());
						if(mapSp2Idx.end() == it1)
							std::cerr << "Error : reaction " << mu << " contains invalid reactant species id '" << speciesReference1->getSpecies() << "'" << std::endl;
						
						SpeciesReference *speciesReference2 = reactionIt.getNextReactant();
						
						// Find the index of the second reactant
						MAP_SP2IDX::iterator	it2 = mapSp2Idx.find(speciesReference2->getSpecies());
						if(mapSp2Idx.end() == it2)
							std::cerr << "Error : reaction " << mu << " contains invalid reactant species id '" << speciesReference2->getSpecies() << "'" << std::endl;
						
						
						// Push the reaction
						
						if(it1->second == it2->second)
						{
							cs.index = it1->second;
							cs.coefficient = -2;
							ptrData->arU.push_back(mu_actual, cs);
							ptrData->arUm.push_back(mu_actual,cs);
							cs.coefficient =  2;
							ptrData->arU.push_back(mu_actual + 1, cs);
							ptrData->arUp.push_back(mu_actual + 1,cs);
						}
						else
						{
							cs.index = it1->second;
							cs.coefficient = -1;
							ptrData->arU.push_back(mu_actual, cs);
							ptrData->arUm.push_back(mu_actual,cs);
							cs.coefficient =  1;
							ptrData->arU.push_back(mu_actual + 1, cs);
							ptrData->arUp.push_back(mu_actual + 1,cs);

							cs.index = it2->second;
							cs.coefficient = -1;
							ptrData->arU.push_back(mu_actual, cs);
							ptrData->arUm.push_back(mu_actual,cs);
							cs.coefficient =  1;
							ptrData->arU.push_back(mu_actual + 1, cs);
							ptrData->arUp.push_back(mu_actual + 1,cs);
						}
					  
						// Add an intermediate species
						cs.index = unIdxAdditionalSpecies;
						cs.coefficient = -1;
						ptrData->arU.push_back(mu_actual, cs);
						ptrData->arUm.push_back(mu_actual,cs);
						cs.coefficient =  1;
						ptrData->arU.push_back(mu_actual + 1, cs);
						ptrData->arUp.push_back(mu_actual + 1,cs);
					  
						mu_actual += 2;
						speciesReference1 = reactionIt.getNextReactant();

						// Find the index of the next reactant
						it1 = mapSp2Idx.find(speciesReference1->getSpecies());
						if(mapSp2Idx.end() == it1)
							std::cerr << "Error : reaction " << mu << " contains invalid reactant species id '" << speciesReference1->getSpecies() << "'" << std::endl;
					  
						// while you do not reach the last reactant
						// create reaction X(i) + Z(n) -> Z(n+1)
						while (reactionIt.hasNext())
						{
							// Add current species to stoichiometry matrices
							cs.index = it1->second;
							cs.coefficient = -1;
							ptrData->arU.push_back(mu_actual, cs);
							ptrData->arUm.push_back(mu_actual,cs);
							cs.coefficient =  1;
							ptrData->arU.push_back(mu_actual + 1, cs);
							ptrData->arUp.push_back(mu_actual + 1,cs);

							// Add an intermediate species
							cs.index = unIdxAdditionalSpecies;
							cs.coefficient = -1;
							ptrData->arU.push_back(mu_actual, cs);
							ptrData->arUm.push_back(mu_actual,cs);
							cs.coefficient =  1;
							ptrData->arU.push_back(mu_actual + 1, cs);
							ptrData->arUp.push_back(mu_actual + 1,cs);

							// Add an intermediate species
							cs.index = unIdxAdditionalSpecies;
							cs.coefficient =  1;
							ptrData->arU.push_back(mu_actual, cs);
							ptrData->arUp.push_back(mu_actual,cs);
							cs.coefficient = -1;
							ptrData->arU.push_back(mu_actual + 1, cs);
							ptrData->arUm.push_back(mu_actual + 1,cs);

							// FIXME Set reaction rate and class(regular)
							for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; si++)
							{
								pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[si];
								sv.ardC[mu_actual]			  = C;
								sv.ardC[mu_actual + 1]		  = (REAL)BREAKUP_REACTION_K;
								sv.arReactClass[mu_actual]	  = pssalib::datamodel::DataModel::D0;
								sv.arReactClass[mu_actual + 1] = pssalib::datamodel::DataModel::D0;
							}

							// Increment reaction counter
							mu_actual+=2;
							unIdxAdditionalSpecies++;
							
							SpeciesReference *speciesReference = reactionIt.getNextReactant();
						
							// Find the index of the second reactant
							MAP_SP2IDX::iterator	it2 = mapSp2Idx.find(speciesReference->getSpecies());
							if(mapSp2Idx.end() == it2)
								std::cerr << "Error : reaction " << mu << " contains invalid reactant species id '" << speciesReference->getSpecies() << "'" << std::endl;
						}
						
						// Add current species to stoichiometry matrices
						cs.index = it1->second;
						cs.coefficient = -1;
						ptrData->arU.push_back(mu_actual, cs);
						ptrData->arUm.push_back(mu_actual,cs);

						// Add an intermediate species
						cs.index = unIdxAdditionalSpecies;
						cs.coefficient = -1;
						ptrData->arU.push_back(mu_actual, cs);
						ptrData->arUm.push_back(mu_actual,cs);
						
					}//breakup
					else
					{
						for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; si++)
						{
						  pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[si];
						  // OVERWRITE THE RATE CONSTANT IF A FILE GIVING THE SPATIAL RATE CONSTANT EXISTS
						  
						  if(ri.spatial_c != NULL)
							{
							  if (ptrSimInfo->bVerbose && si == 0)
								std::cout << "Reaction " << mu << " rate constant is spatial";
							  sv.ardC[mu_actual] = ri.spatial_c[si];
							}

						  //END OVERWRITE 
						  else{
							sv.ardC[mu_actual] = (bReverseLoop ? ri.c_r : ri.c);
						  }
							sv.arReactClass[mu_actual] = pssalib::datamodel::DataModel::D0;
						}

						// Fill in the stoichiometry matrix
						for (UNSIGNED_INTEGER j = 0; j < unEntries; j++)
						{
							SpeciesReference *speciesReference = bReverseLoop ? reaction->getProduct(j) : reaction->getReactant(j);
							
							if (ptrSimInfo->bVerbose)
								std::cout << "reactant " << j << " : species id = '" << speciesReference->getSpecies() << "' : ";

							MAP_SP2IDX::iterator	it = mapSp2Idx.find(speciesReference->getSpecies());
							if(mapSp2Idx.end() == it)
								std::cerr << "Error : reaction " << mu_actual << " contains invalid reactant species id '" << speciesReference->getSpecies() << "'" << std::endl;

							cs.index = it->second;
							cs.coefficient = -1 * ((INTEGER)fabs(speciesReference->getStoichiometry()));
							ptrData->arU.push_back(mu_actual, cs);
							ptrData->arUm.push_back(mu_actual,cs);
						}
						
						if (ptrSimInfo->bVerbose)
							std::cout << std::endl;
					}
				}

				unEntries = bReverseLoop ? reaction->getNumReactants() : reaction->getNumProducts();
				if(unEntries == 0)
				{
					std::cerr << "Warning : reaction " << reaction->getId() << (bReverseLoop ? "(reverse)" : "") << " has no products! Assuming 'reservoir'." << std::endl;
					cs.index = 0;
					cs.coefficient = 1;
					ptrData->arU.push_back(mu_actual, cs);
					ptrData->arUp.push_back(mu_actual,cs);
				}
				else
				{
					if (ptrSimInfo->bVerbose)
						std::cout << "Reaction " << reaction->getId() << " has " << reaction->getNumProducts() << " products : " << std::endl << '\t';

					for (UNSIGNED_INTEGER j = 0; j < unEntries; ++j)
					{
						SpeciesReference *speciesReference = bReverseLoop ? reaction->getReactant(j) : reaction->getProduct(j);

						if (ptrSimInfo->bVerbose)
							std::cout << "product " << j << " : species id = '" << speciesReference->getSpecies() << "' : ";

						MAP_SP2IDX::iterator it = mapSp2Idx.find(speciesReference->getSpecies());
						if(mapSp2Idx.end() == it)
							std::cerr << "reaction " << mu_actual << " contains invalid product species id '" << speciesReference->getSpecies() << "'" << std::endl;

						cs.index = it->second;
						cs.coefficient = ((INTEGER)fabs(speciesReference->getStoichiometry()));
						ptrData->arU.push_back(mu_actual, cs);
						ptrData->arUp.push_back(mu_actual,cs);
					}

					if (ptrSimInfo->bVerbose)
						std::cout << std::endl;
				}

				//
				if(!bReverseLoop&&reaction->getReversible())
				{
					bReverseLoop = true;
					--mu;
				}
				else if(bReverseLoop)
					bReverseLoop = false;
			}

			// Create reduced version of stochimetry matrix
			
			for (int i = 0 ; i < ptrData->arU.get_rows() ; i++)
			{
				for (int j = 0 ; j < ptrData->arU.get_cols(i) ; j++)
				{
					// search if index has been already pushed
				  
					int k;
					for (k = 0 ; k < ptrData->arU_r.get_cols(i) ; k++)
					{
						if (ptrData->arU(i,j).index == ptrData->arU(i,k).index)
						{
							// founded, update coefficient
							
							ptrData->arU_r(i,k).coefficient += ptrData->arU(i,j).coefficient;
						  
							break;
						}
					}
					
					// push if not finded
					
					if (k == ptrData->arU_r.get_cols(i))
						ptrData->arU_r.push_back(i,ptrData->arU(i,j));
				}
			}
			
			// Create reduced version of stochimetry matrix
			
			for (int i = 0 ; i < ptrData->arUm.get_rows() ; i++)
			{
				for (int j = 0 ; j < ptrData->arUm.get_cols(i) ; j++)
				{
					// search if index has been already pushed
				  
					int k;
					for (k = 0 ; k < ptrData->arUm_r.get_cols(i) ; k++)
					{
						if (ptrData->arUm(i,j).index == ptrData->arUm(i,k).index)
						{
							// founded, update coefficient
							
							ptrData->arUm_r(i,k).coefficient += ptrData->arUm(i,j).coefficient;
						  
							break;
						}
					}
					
					// push if not finded
					
					if (k == ptrData->arUm_r.get_cols(i))
						ptrData->arUm_r.push_back(i,ptrData->arUm(i,j));
				}
			}
			
			// Create reduced version of stochimetry matrix
			
			for (int i = 0 ; i < ptrData->arUp.get_rows() ; i++)
			{
				for (int j = 0 ; j < ptrData->arUp.get_cols(i) ; j++)
				{
					// search if index has been already pushed
				  
					int k;
					for (k = 0 ; k < ptrData->arUp_r.get_cols(i) ; k++)
					{
						if (ptrData->arUp(i,j).index == ptrData->arUp(i,k).index)
						{
							// founded, update coefficient
							
							ptrData->arUp_r(i,k).coefficient += ptrData->arUp(i,j).coefficient;
						  
							break;
						}
					}
					
					// push if not finded
					
					if (k == ptrData->arUp_r.get_cols(i))
						ptrData->arUp_r.push_back(i,ptrData->arUp(i,j));
				}
			}
			
			return true;
		}

		double GroupingModule::getDiffusion(XMLNode * xn, std::string str)
		{
			double diff;
		  
			if ( xn->isEOF() )
			{
				// root node is a dummy node
				for ( int i = 0; i < xn->getNumChildren(); i++ )
				{
					// access to each child node of the dummy node.
					XMLNode& xnChild = xn->getChild(i);
					
					if (xnChild.getName() == str)
					{
						const XMLAttributes atr = xnChild.getAttributes();
						std::istringstream i(xnChild.getAttrValue(atr.getIndex("diffusion")));
						
						i >> diff;
						
						return diff;
					}
				}
			}                                                                                      
			else                                                                                   
			{                                                                                      
				// root node is NOT a dummy node
				
				std::string tmp = xn->getName();
				
				if (tmp.compare(std::string("annotation")) != 0 || xn->getNumChildren() == 0)
					return 0.0;
				
				XMLNode& xnChild = xn->getChild(0);
				
				if (xnChild.getName().compare(std::string("diffusion")) != 0 || xn->getNumChildren() == 0)
					return 0.0;
				
				for ( int i = 0; i < xnChild.getNumChildren(); i++ )
				{
					// access to each child node
					XMLNode& xnCR = xnChild.getChild(i);
					
					if (xnCR.getName() == str)
					{
						const XMLAttributes atr = xnCR.getAttributes();
						std::istringstream i(xnCR.getAttrValue(atr.getIndex("diffusion")));
						
						i >> diff;
						
						return diff;
					}
				}
			}
		}
		
		bool GroupingModule::fill_diffusion(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Get the model
			Model *model = ptrSimInfo->getSBMLModel();
			pssalib::datamodel::DataModel * ptrData = ptrSimInfo->getDataModel();

			XMLNode * nd = model->getListOfSpecies()->getAnnotation();
		  
			if ( nd == NULL && ptrData->unSubvolumes != 1)
			{
				std::cerr << "Error : diffusion constant has to be defined into annotation tag \n" << std::endl;
				return false;
			}
			
			if (ptrData->unSubvolumes != 1)
			{
				for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes ; si++)
				{
					for(UNSIGNED_INTEGER nu = 0; nu < ptrData->unSpecies - 1; nu++)
					{

						ptrData->arSubvolumes[si].ardD[nu] = getDiffusion(nd,std::string("d") + model->getSpecies(nu)->getName());
					}
				}
			}
			else
			{
				// One volume no diffusion
			  
				for(UNSIGNED_INTEGER nu = 0; nu < ptrData->unSpecies - 1; nu++)
				{
					ptrData->arSubvolumes[0].ardD[nu] = 0.0;
				}
			}
			  

			return true;
		}

		int factorial(int f)
		{
			int a = 1;
		  
			for (int i = f ; i > 0 ; i--)
			{
				a *= i;
			}
			
			return a;
		}
		
		bool GroupingModule::adjust_reaction_rates_by_volume(pssalib::datamodel::DataModel* ptrData)
		{
			for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; si++)
			{
				pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[si];
				for(UNSIGNED_INTEGER mu = 0; mu < ptrData->unReactions; mu++)
				{
					UNSIGNED_INTEGER num_reactants = 0;
					bool same_species = false;

					for(UNSIGNED_INTEGER nu = 0; nu < ptrData->arUm.get_cols(mu); nu++)
					{
						pssalib::datamodel::DataModel::ChemicalSpecies& cs = ptrData->arUm(mu, nu);
						if (cs.index != 0)
						{
							num_reactants += -ptrData->arUm(mu, nu).coefficient;
						}
						if (cs.coefficient == -2)
						{
							same_species = true;
						}
					}

					REAL factor;
					
					factor = pow(sv.dOmega,(int)(1-num_reactants));
					
					//    __
					//    || !stochimetry(nu)
					//   nu = species involved
					//
					
					for(UNSIGNED_INTEGER nu = 0; nu < ptrData->arUm.get_cols(mu); nu++)
					{
						factor *= factorial(abs(ptrData->arUm(mu,nu).coefficient));
					}
/*					switch (num_reactants)
					{
					case 0:
						factor = sv.dOmega;
						break;

					case 1:
						factor = 1.0;
						break;

					case 2:
						factor = 1.0 / sv.dOmega;
						if (same_species)
							factor *= 2.0;
						break;

					default:
						;
					}*/

					sv.ardC[mu] *= factor;
				}

				for(UNSIGNED_INTEGER nu = 0; nu < ptrData->unSpecies-1; nu++)
				{
					sv.ardD[nu] /= sv.dH * sv.dH;
				}
			}

			return true;
		}

	  //////////////////////////////:
	  double* GroupingModule::getSpatialC(XMLNode * annotation,int gridpoints)
	  {
		int n_subvolumes = gridpoints * gridpoints;
		  
		XMLNode& child = annotation->getChild(0);
		const XMLAttributes atr = child.getAttributes();
		std::string content = child.getAttrValue(atr.getIndex("rate"));
		
		
		// Get the content of the xml node and split it...
		std::vector<std::string> strs;  
		boost::split(strs, content, boost::is_any_of(";"));
		
		
		if(strs.size()!= n_subvolumes)
		  {
			std::cout << "BAD NUMBER OF COEFFICIENTS"<< strs.size() << "!="<< n_subvolumes << std::endl;
			std::cout << content << std::endl;
		  }
		  
		// Cast the elements into an array of double...
		double *spatial_c;
		spatial_c = new double[n_subvolumes];
		for (int i = 0; i < n_subvolumes; ++i)
		  {
			spatial_c[i] =   boost::lexical_cast<double>(strs[i]); 
		  }

		return spatial_c;
	  }

	  ///////////////////////////////////////////

	  
	}
}


