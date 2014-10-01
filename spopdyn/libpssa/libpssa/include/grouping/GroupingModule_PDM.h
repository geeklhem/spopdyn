/**
 * \file GroupingModule_PDM.h
 *
 * \copydoc GroupingModule.h
 *
 * \date   13.12.2010
 * \author sanyi
 */

#ifndef GROUPINGMODULE_PDM_H_
#define GROUPINGMODULE_PDM_H_

#include "./GroupingModule.h"

namespace pssalib
{
	namespace grouping
	{
		class GroupingModule_PDM : public GroupingModule
		{
		////////////////////////////////
		// Constructors
		public:
			// Default constructor
			GroupingModule_PDM();

			// Copy constructor
			GroupingModule_PDM(GroupingModule &);

			// Destructor
	virtual ~GroupingModule_PDM();

		////////////////////////////////
		// Methods
		protected:
			// Check if break-up of non-elementary reactions is required
	virtual	bool	do_breakup() const;

		public:
			// Initialize data structures (called before each trial)
	virtual bool initialize(pssalib::datamodel::DataModel* ptrData);
		};
	}
}

#endif /* GROUPINGMODULE_PDM_H_ */
