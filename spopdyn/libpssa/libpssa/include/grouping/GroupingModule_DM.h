/**
 * \file GroupingModule_DM.h
 *
 * \copydoc GroupingModule.h
 *
 * \date   13.12.2010
 * \author sanyi
 */

#ifndef GROUPINGMODULE_DM_H_
#define GROUPINGMODULE_DM_H_

#include "./GroupingModule.h"

namespace pssalib
{
	namespace grouping
	{
		class GroupingModule_DM : public GroupingModule
		{
		////////////////////////////////
		// Constructors
		public:
			// Default constructor
			GroupingModule_DM();

			// Copy constructor
			GroupingModule_DM(GroupingModule &);

			// Destructor
	virtual ~GroupingModule_DM();

		////////////////////////////////
		// Methods
		protected:
			// Check if break-up of non-elementary reactions is required
	virtual	bool	do_breakup() const;

		public:
			// Calculate total propensity
	virtual bool initialize(pssalib::datamodel::DataModel* ptrData);
		};
	}
}

#endif /* GROUPINGMODULE_DM_H_ */
