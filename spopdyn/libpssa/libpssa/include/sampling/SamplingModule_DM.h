/**
 * \file SamplingModule_DM.h
 * \copybrief SamplingModule.h
 *
 * \date   13.12.2010
 * \author sanyi
 */

#ifndef SAMPLINGMODULE_DM_H_
#define SAMPLINGMODULE_DM_H_

#include "./SamplingModule.h"

namespace pssalib
{
	namespace sampling
	{
		/**
		 * \class SamplingModule_DM
		 * \brief Purpose of this class is to sample reaction index and time for the Gillespie's Direct Method.
		 *
		 * \copydetails SamplingModule
		 */
		class SamplingModule_DM : public SamplingModule
		{
		/////////////////////////////////////
		// Constructors
		public:
			// Default Constructor
			SamplingModule_DM();
			// Destructor
	virtual	~SamplingModule_DM();

		//////////////////////////////
		// Methods
		protected:
			// Sample next reaction index
	virtual bool sample_reaction(pssalib::datamodel::SimulationInfo* ptrSimInfo);
		};
	}
}

#endif /* GROUPINGMODULE_DM_H_ */
