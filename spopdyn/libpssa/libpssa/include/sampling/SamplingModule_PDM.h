/**
 * \file SamplingModule_PDM.h
 * \copybrief SamplingModule.h
 *
 * \date   13.12.2010
 * \author sanyi
 */

#ifndef SAMPLINGMODULE_PDM_H_
#define SAMPLINGMODULE_PDM_H_

#include "./SamplingModule.h"

namespace pssalib
{
	namespace sampling
	{
		/**
		 * \class SamplingModule_PDM
		 * \brief Purpose of this class is to sample reaction index and time for the Partial Propensity Direct Method.
		 *
		 * \copydetails SamplingModule
		 */
		class SamplingModule_PDM : public SamplingModule
		{
		/////////////////////////////////
		// Constructors
		public:
			// Constructor
			SamplingModule_PDM();
			// Destructor
	virtual	~SamplingModule_PDM();

		////////////////////////////////
		// Sampling module methods
		protected:
			// Sample reaction index
	virtual bool sample_reaction(pssalib::datamodel::SimulationInfo* ptrSimInfo);
		};
	}	// namespace sampling
}	// namespace pssalib

#endif /* SAMPLINGMODULE_PDM_H_ */
