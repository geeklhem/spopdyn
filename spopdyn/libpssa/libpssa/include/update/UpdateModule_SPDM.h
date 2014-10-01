/**
 * \file UpdateModule_SPDM.h
 * \copybrief UpdateModule.h
 *
 * \date   14.15.2011
 * \author sanyi
 */

#ifndef UPDATEMODULE_SPDM_H_
#define UPDATEMODULE_SPDM_H_

#include "./UpdateModule_PDM.h"

namespace pssalib
{
	namespace update
	{
		/**
		 * \class UpdateModule_SPDM
		 * \brief Purpose of this class is to update the datastructures for the Partial Propensity Direct Method after a reaction has fired.
		 *
		 * \copydetails UpdateModule
		 */
		class UpdateModule_SPDM : public UpdateModule_PDM
		{
		////////////////////////////////
		// Constructors
		public:
			// Default Constructor
			UpdateModule_SPDM();
			// Destructor
	virtual	~UpdateModule_SPDM();

		////////////////////////////////
		// Update module methods
		protected:
			// Update structures
	virtual bool update_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo);
		};
	}	// namespace update
}	// namespace pssalib

#endif /* UPDATEMODULE_SPDM_H_ */
