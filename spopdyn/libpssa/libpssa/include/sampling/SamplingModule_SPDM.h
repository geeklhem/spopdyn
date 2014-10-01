/**
 * \file SamplingModule_SPDM.h
 * \copybrief SamplingModule.h
 *
 * \date   14.05.2011
 * \author sanyi
 */

#ifndef SAMPLINGMODULE_SPDM_H_
#define SAMPLINGMODULE_SPDM_H_

#include "./SamplingModule_PDM.h"

namespace pssalib
{
	namespace sampling
	{
		/**
		 * \class SamplingModule_SPDM
		 * \brief Purpose of this class is to sample reaction index and time for the Sorting Partial Propensity Direct Method.
		 *
		 * \copydetails SamplingModule
		 */
		class SamplingModule_SPDM : public SamplingModule_PDM
		{
		/////////////////////////////////
		// Constructors
		public:
			// Constructor
			SamplingModule_SPDM();
			// Destructor
	virtual	~SamplingModule_SPDM();

		////////////////////////////////
		// Sampling module methods
		protected:
			// Sample reaction index
	virtual bool sample_reaction(pssalib::datamodel::SimulationInfo* ptrSimInfo);
		};
	}	// namespace sampling
}	// namespace pssalib

#endif /* SAMPLINGMODULE_SPDM_H_ */
