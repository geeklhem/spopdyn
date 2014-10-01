/**
 * \file UpdateModule_DM.h
 * \copybrief UpdateModule.h
 *
 * \date   15.12.2010
 * \author sanyi
 */


#ifndef UPDATEMODULE_DM_H_
#define UPDATEMODULE_DM_H_

#include "./UpdateModule.h"
#include "../../include/datamodel/DataModel_DM.h"

namespace pssalib
{
	namespace update
	{
		/**
		 * \class UpdateModule_DM
		 * \brief Purpose of this class is to update the datastructures for the Gillespie's Direct Method after a reaction has fired.
		 *
		 * \copydetails UpdateModule
		 */
		class UpdateModule_DM : public UpdateModule
		{
		////////////////////////////////
		// Constructors
		public:
		  
			// Constructor
			UpdateModule_DM();
			// Destructor
	virtual ~UpdateModule_DM();

		////////////////////////////////
		// Update module methods
		protected:
			// Update structures
	virtual bool update_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo);
	
		private:
	  
	void update_structures_species(pssalib::datamodel::DataModel_DM* ptrDMData, UNSIGNED_INTEGER svix, UNSIGNED_INTEGER species, INTEGER species_coeff, REAL& dlt_a_0);
	
		};
	  
	}
}

#endif /* UPDATEMODULE_DM_H_ */
