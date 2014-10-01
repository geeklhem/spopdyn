/**
 * \file UpdateModule_PSSACR.h
 * \copybrief UpdateModule.h
 *
 * \date   22.01.2011
 * \author sanyi
 */

#ifndef UPDATEMODULE_PSSACR_H_
#define UPDATEMODULE_PSSACR_H_

#include "./UpdateModule_PDM.h"

namespace pssalib
{
	namespace datamodel
	{
		class DataModel_PSSACR;
	}

	namespace update
	{
		/**
		 * \class UpdateModule_PSSACR
		 * \brief Purpose of this class is to update the datastructures for the Partial Propensity Direct Method with Composition-Rejection sampling after a reaction has fired.
		 *
		 * \copydetails UpdateModule
		 */
		class UpdateModule_PSSACR : public UpdateModule_PDM
		{
		////////////////////////////////
		// Update module attributes
		protected:
			//! Stores indexes of all changed Sigmas
			std::vector<UNSIGNED_INTEGER>	vIdxSigmaChanged;

		////////////////////////////////
		// Constructors
		public:
			// Constructor
			UpdateModule_PSSACR();
			// Destructor
	virtual	~UpdateModule_PSSACR();

		////////////////////////////////
		// Update module methods
		protected:
			// Update structures
	virtual bool update_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo);

		private:
			void update_structures_species(pssalib::datamodel::DataModel_PSSACR* ptrPSRDCRData, UNSIGNED_INTEGER svix, UNSIGNED_INTEGER species);
		};
	}	// namespace update
}	// namespace psrdlib

#endif /* UPDATEMODULE_PSSACR_H_ */
