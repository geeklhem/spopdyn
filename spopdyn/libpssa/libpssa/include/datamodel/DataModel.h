/**
 * \file DataModel.h
 * \brief Declares a datatype containing all the datastructures commonly used in all SSAs
 *
 * \details \c DataModel class stores the stoichiometry of the reaction network,
 * amounts of species, reaction parameters (propensities, classes).
 * It also contains some datatypes commonly used by the library internally, e.g.,
 * \c ChemicalSpecies, \c ReactionClass, \c DelayedReaction.
 * This class does \b not contain any methods to operate on the data.
 * It is used as a base class for classes implementing datastructures for
 * different SSA implementations included in the library.
 *
 * \date   14.12.2010
 * \author sanyi
 */

#ifndef DATAMODEL_H_
#define DATAMODEL_H_

#include "../typedefs.h"
#include "CompositionRejectionSamplerData.h"

namespace pssalib
{
	namespace datamodel
	{
		////////////////////////////////////////
		//! Generic data model
		class DataModel
		{
		public:
		  
			std::vector<std::string> speciesName;
		  
			//! Data structure for reactants and products
			typedef struct tagChemicalSpecies
			{
				UNSIGNED_INTEGER	index;		//!< index of species
				INTEGER				coefficient;//!< coefficient with which it appears in reaction

				tagChemicalSpecies() : index(0), coefficient(0)	{	};
				~tagChemicalSpecies()	{	};

		inline	bool	operator<(const tagChemicalSpecies& right) const
				{
					return (coefficient < right.coefficient);
				};

		inline	bool	operator==(const tagChemicalSpecies& right) const
				{
					return (coefficient == right.coefficient);
				};

				friend std::ostream & operator<<( std::ostream & output, const tagChemicalSpecies& C)
				{
					output << "[I = " << C.index << ", K = " << C.coefficient << "]";
					return output;
				};
			} ChemicalSpecies;

			/*!
			 * \enum psrdlib::datamodel::DataModel::ReactionClass
			 * \brief Enum for reactions classes. The abbreviations are taken from (Cai, 2007)
			 *
			 *\details  This enum provides three named constants, each representing a reaction class supported by the library.
			 * They are used during initializain the \b Update module to determine the update mode.
			 */
			typedef enum tagReactionClass
			{
				/*! Ordinary chemical reaction. */
				D0,
				/*! Delayed reaction (class 1) : non-consuming reaction.
				 *  \attention FIXME Not tested!
				 */
				D1,
				/*! Delayed reaction (class 2) : consuming reaction. */
				D2
			} ReactionClass;

			//! Struct for delayed reactions.
			//! FIXME Check if this fits both D1 & D2 reaction classes.
			typedef struct tagDelayedReaction
			{
				UNSIGNED_INTEGER	index;		//!< index of reaction
				REAL				time;		//!< time that passes before it fires

				//! Default constructor
				tagDelayedReaction() : index(0), time(0.0)	{	};

				//! Constructor
				tagDelayedReaction(UNSIGNED_INTEGER i, REAL t) : index(i), time(t)	{	};

				//! Destructor
				~tagDelayedReaction()	{	};

		inline	bool	operator<(const tagDelayedReaction& right) const
				{
					return (time < right.time);
				};
			} DelayedReaction;

			struct Subvolume
			{
				Subvolume();
				~Subvolume();

				UNSIGNED_INTEGER				volumeIndex;

				// Species
				//

				//! Vector of current species population
				INTEGER							*aruN;

				//! Total propensity
				REAL							dTotalPropensity;

				// Reactions
				//

				//! Number of diffusion reactions
				UNSIGNED_INTEGER				unDiffusionReactions;

				//! 
				UNSIGNED_INTEGER				*arNeighborVolumes;

				//! Rate constants
				REAL							*ardC;

				//! Diffusion constants
				REAL							*ardD;

				//! Vector of reaction classes
				ReactionClass					*arReactClass;

				//! Reaction delays
				REAL 							*arDelay;

				// Subvolume description
				REAL							dH;
				REAL							dOmega;
			};
			
			class SubvolumesInit
			{
				private:
					
					int unSpecies;
					int nsub;
					Subvolume * subVol;
					
				public:
				  
					SubvolumesInit(int nsp, int nsb, Subvolume * subVolt)	{unSpecies = nsp; nsub = nsb; subVol = subVolt;};
					~SubvolumesInit()	{};
				  
					int getNumSub()
					{
						return nsub;
					}
					
					void setInit(int v, int sp ,int init_am)
					{
						if (v >= nsub)
							throw "SubvolumesInit: Subvolumes OutofBound";
					  
						if (sp >= unSpecies)
							throw "SubvolumesInit: OutofBound species";
						
						subVol[v].aruN[sp] = init_am;
					}
			};

		public:
			DataModel();
			DataModel(DataModel &);
	virtual	~DataModel();
	virtual void Cleanup();
	inline  bool IsDiffusionReaction(UNSIGNED_INTEGER ix) { return ix >= unReactions; }

		public:
			//! \internal Model name
			STRING							strModelName;

			// Stoichiometry
			//

			//! Stoichiometry matrix representation (matrices U1 and U2 from the paper)
			//! arU	- stoichiometry info for all reactions
			MatrixVarLen<ChemicalSpecies>	arU; 
			//! Stoichiometry matrix representation (matrices U1 and U2 from the paper)
			//! arUp	- stoichiometry info of products for all reactions
			MatrixVarLen<ChemicalSpecies>	arUp;
			//! Stoichiometry matrix representation (matrices U1 and U2 from the paper)
			//! arUm	- all stoichiometry info of reactants for all reactions
			MatrixVarLen<ChemicalSpecies>	arUm;
			
			//! Stoichiometry matrix representation (matrices U1 and U2 from the paper) reduced version
			//! arU	- stoichiometry info for all reactions
			MatrixVarLen<ChemicalSpecies>	arU_r; 
			//! Stoichiometry matrix representation (matrices U1 and U2 from the paper) reduced version
			//! arUp	- stoichiometry info of products for all reactions
			MatrixVarLen<ChemicalSpecies>	arUp_r;
			//! Stoichiometry matrix representation (matrices U1 and U2 from the paper) reduced version
			//! arUm	- all stoichiometry info of reactants for all reactions
			MatrixVarLen<ChemicalSpecies>	arUm_r;
			
			UNSIGNED_INTEGER				unSpecies;
			UNSIGNED_INTEGER				unReactions;
			UNSIGNED_INTEGER				unSubvolumes;
			Subvolume						*arSubvolumes;
			UNSIGNED_INTEGER				unMiddleVolumeIx;

			CompositionRejectionSamplerData crsdVolume;

			//! Total propensity
			REAL							dTotalPropensity;

		/////////////////////////////////////
		// Sampling variables
		public:
			//! Sampled reaction index
			INTEGER							mu;

			//! Sampled volume
			UNSIGNED_INTEGER				sv;

			//! Sampled destination volume (for diffusion reactions)
			UNSIGNED_INTEGER				svDst;

			//! Queued reactions (both D1 & D2 class reactions)
			std::vector<DelayedReaction>	vQueuedReactions;
		};
	}
}

#endif /* DATAMODEL_H_ */
