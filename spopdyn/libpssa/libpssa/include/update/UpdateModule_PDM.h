/**
 * \file UpdateModule_PDM.h
 * \copybrief UpdateModule.h
 *
 * \date   15.12.2010
 * \author sanyi
 */

#ifndef UPDATEMODULE_PDM_H_
#define UPDATEMODULE_PDM_H_

#include "./UpdateModule.h"

namespace pssalib
{
	namespace datamodel
	{
		class DataModel_PDM;
	}

	namespace update
	{
		/**
		 * \class UpdateModule_PDM
		 * \brief Purpose of this class is to update the datastructures for the Partial Propensity Direct Method after a reaction has fired.
		 *
		 * \copydetails UpdateModule
		 */
		class UpdateModule_PDM : public UpdateModule
		{
		////////////////////////////////
		// Constructors
		public:
			// Default Constructor
			UpdateModule_PDM();
			// Destructor
	virtual	~UpdateModule_PDM();
	
			static double comb(int n, int coeff);

		////////////////////////////////
		// Update module methods
		protected:
			// Update structures
	virtual bool update_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo);

		private:
			void update_structures_species(pssalib::datamodel::DataModel_PDM* ptrPDMData, UNSIGNED_INTEGER svix, UNSIGNED_INTEGER species, INTEGER species_coeff, REAL& dlt_a_0);
		};
	}	// namespace update
}	// namespace pssalib

#endif /* UPDATEMODULE_PDM_H_ */
