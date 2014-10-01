/**
 * \file GroupingModule_SPDM.h
 *
 * \copydoc GroupingModule.h
 *
 * \date   13.12.2010
 * \author sanyi
 */

#ifndef GROUPINGMODULE_SPDM_H_
#define GROUPINGMODULE_SPDM_H_

#include "./GroupingModule_PDM.h"

namespace pssalib
{
	namespace grouping
	{
		class GroupingModule_SPDM : public GroupingModule_PDM
		{
		////////////////////////////////
		// Constructors
		public:
			// Default constructor
			GroupingModule_SPDM();

			// Copy constructor
			GroupingModule_SPDM(GroupingModule &);

			// Destructor
	virtual ~GroupingModule_SPDM();

		////////////////////////////////
		// Methods
		public:
			// Initialize data structures (called before each trial)
	virtual bool initialize(pssalib::datamodel::DataModel* ptrData);
		};
	}
}

#endif /* GROUPINGMODULE_SPDM_H_ */
