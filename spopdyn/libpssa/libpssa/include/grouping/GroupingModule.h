/**
 * \file GroupingModule.h
 * \brief Grouping module is responsible for loading the data from an SBML model and filling in the data structures.
 *
 * \c GroupingModule class internally stores a mapping from species identifier
 * to species index used in all of the library's datastructures.
 * It is used as a base class for classes implementing grouping for
 * different SSA implementations included in the library.
 *
 * \date   13.12.2010
 * \author sanyi
 */

#ifndef GROUPINGMODULE_H_
#define GROUPINGMODULE_H_

#include "../typedefs.h"
#include "../datamodel/DataModel.h"
#include "../datamodel/SimulationInfo.h"

// Forward declarations
namespace pssalib
{
	namespace datamodel
	{
		class DataModel;
		class SimulationInfo;
	}
}

namespace pssalib
{
	namespace grouping
	{
		class GroupingModule
		{
		private:
		
			bool eval_error;
		  
			double eval(const ASTNode * nd, double v1, double v2);
			double evaluateMath(const ASTNode * Mnode, KineticLaw * L, Model * model);
			double searchPar(std::string str, KineticLaw * L, Model * model);
			
		protected:
		  
			//! Temporary storage class (SBML=>PSSAlib translation)
			typedef struct tagReactionInfo
			{
				pssalib::datamodel::DataModel::ReactionClass	rc;
				//! Specific probability rate (forward & reverse) and delay
				REAL				c, c_r, d;
				// Number of additional reactions (forward & reverse)
			  UNSIGNED_INTEGER	n, n_r;

			  // Spatial rate
			  double* spatial_c;
			  //
			  
				tagReactionInfo() : rc(pssalib::datamodel::DataModel::D0),
						c(0.0), c_r(0.0), d(0.0), n(0), n_r(0)
				{

				};

				void clear()
				{
					rc	= pssalib::datamodel::DataModel::D0;
					c	= 0.0;
					c_r	= 0.0;
					d	= 0.0;
					n	= 0;
					n_r	= 0;
				};
			} ReactionInfo;

		////////////////////////////////
		// Constructors
		public:
			// Default constructor
			GroupingModule();

			// Copy constructor
			GroupingModule(GroupingModule &);

			// Destructor
	virtual ~GroupingModule();

		////////////////////////////////
		// Methods
		protected:
			// Check if break-up of non-elementary reactions is required
	virtual	bool do_breakup() const = 0;

		public:
			// Read the SBML model
	virtual bool preinitialize(pssalib::datamodel::SimulationInfo* ptrSimInfo);

			// Initialize data structures (called before each trial)
	virtual bool initialize(pssalib::datamodel::DataModel* ptrData);

			// Initialize composition-rejection samples for subvolumes.
	virtual bool postinitialize(pssalib::datamodel::DataModel* ptrData);

		private:
			bool parse_parameters(pssalib::datamodel::SimulationInfo* ptrSimInfo, Model* model, std::vector<ReactionInfo>& vReactionInfo, REAL volume);
			bool setup_volumes_single(pssalib::datamodel::DataModel* ptrData, pssalib::datamodel::SimulationInfo::BoundaryConditionsType bc, REAL omega);
			bool setup_volumes_homogeneous_2d(pssalib::datamodel::DataModel* ptrData, pssalib::datamodel::SimulationInfo::BoundaryConditionsType bc, UNSIGNED_INTEGER numGridPoints, REAL omega);
			bool setup_volumes_homogeneous_3d(pssalib::datamodel::DataModel* ptrData, pssalib::datamodel::SimulationInfo::BoundaryConditionsType bc, UNSIGNED_INTEGER numGridPoints, REAL omega);
			bool initialize_volumes(pssalib::datamodel::DataModel* ptrData, UNSIGNED_INTEGER numDiffusionReactions, REAL omega, REAL h);
			bool initialize_mappings(pssalib::datamodel::SimulationInfo* ptrSimInfo, pssalib::datamodel::DataModel* ptrData, Model* model);
			bool fill_mappings(pssalib::datamodel::SimulationInfo* ptrSimInfo, pssalib::datamodel::DataModel* ptrData, Model* model, const std::vector<ReactionInfo>& vReactionInfo);
			bool fill_diffusion(pssalib::datamodel::SimulationInfo* ptrSimInfo);
			bool adjust_reaction_rates_by_volume(pssalib::datamodel::DataModel* ptrData);
			void printsp(int coef, std::string name);
			void PrintSimulationInformation(pssalib::datamodel::SimulationInfo* ptrSimInfo);
			double getDiffusion(XMLNode * xn, std::string str);
			double* getSpatialC(XMLNode * annotation, int gridpoints);
			
		public:
			//! Define a type for index=>sbml_species mapping
			typedef boost::unordered_map	<STRING,	UNSIGNED_INTEGER>	MAP_SP2IDX;
			//! Define a type for sbml_species=>index mapping
			typedef	std::vector	<STRING>									VEC_IDX2SP;

			//! sbml_species=>index mapping
			MAP_SP2IDX	mapSp2Idx;

			//! index=>sbml_species mapping
			VEC_IDX2SP	arIdx2Sp;

		protected:

			//! Flag for successful loading of data
			bool	bDataLoaded;

			UNSIGNED_INTEGER unNumAdditionalTotal;
			UNSIGNED_INTEGER unReverseReactions;
		};
	}
}

#endif /* GROUPINGMODULE_H_ */
