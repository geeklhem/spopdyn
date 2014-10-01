/**
 * \file GroupingModule_PSSACR.h
 *
 * \copydoc GroupingModule.h
 *
 * \date   22.02.2011
 * \author sanyi
 */

#ifndef GROUPINGMODULE_PSSACR_H_
#define GROUPINGMODULE_PSSACR_H_

#include "./GroupingModule_PDM.h"

namespace pssalib
{
	namespace grouping
	{
		class GroupingModule_PSSACR : public GroupingModule_PDM
		{
		////////////////////////////////
		// Constructors
		public:
			// Default constructor
			GroupingModule_PSSACR();

			// Copy constructor
			GroupingModule_PSSACR(GroupingModule &);

			// Destructor
	virtual ~GroupingModule_PSSACR();

		////////////////////////////////
		// Methods
		public:
			// Initialize data structures (called before each trial)
	virtual bool initialize(pssalib::datamodel::DataModel* ptrData);
		};
	}
}

#endif /* GROUPINGMODULE_PSSACR_H_ */
