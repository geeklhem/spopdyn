/**
 * \file UpdateModule.h
 * \brief \b Update module is responsible for updating the data structures after a reaction has fired.
 *
 * \date   15.12.2010
 * \author sanyi
 */

#ifndef UPDATEMODULE_H_
#define UPDATEMODULE_H_

#include "../typedefs.h"

// Forward declarations
namespace pssalib
{
	namespace datamodel
	{
		class SimulationInfo;
	}
}

namespace pssalib
{
	namespace update
	{
		/**
		 * \class UpdateModule
		 * \brief This class serves as a base class for all update schemes used within the library.
		 *
		 * \details Provides the functionality necessary to update the respective data structures after a reaction has fired.
		 * All protected methods of this class may be overloaded, offering a flexible way to extend its functionality.
		 */
		class UpdateModule
		{
		/////////////////////////////////
		// Constructors
		public:
			// Constructor
			UpdateModule();
			// Destructor
	virtual ~UpdateModule();

		/////////////////////////////////
		// Update module methods
		public:
			// Perform the update step
	virtual bool do_update(pssalib::datamodel::SimulationInfo* ptrSimInfo);

		protected:
			// Schedule a delayed reaction
	virtual bool schedule_delayed(pssalib::datamodel::SimulationInfo* ptrSimInfo);

			// Update population
	virtual bool update_species(pssalib::datamodel::SimulationInfo* ptrSimInfo);

			//! Update structures
	virtual bool update_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo) = 0;
			bool update_volume_structures(pssalib::datamodel::SimulationInfo* ptrSimInfo);
		};
	}	// namespace update
}	// namespace pssalib

#endif /* UPDATEMODULE_H_ */
