/**
 * \file SamplingModule_PSSACR.h
 * \copybrief SamplingModule.h
 *
 * \date   22.02.2011
 * \author sanyi
 */

#ifndef SAMPLINGMODULE_PSSACR_H_
#define SAMPLINGMODULE_PSSACR_H_

#include "./SamplingModule_PDM.h"
#include "./CompositionRejectionSampler.h"

namespace pssalib
{
	namespace sampling
	{
		/**
		 * \class SamplingModule_PSSACR
		 * \brief Purpose of this class is to sample reaction index and time for the Partial Propensity Direct Method with Compositio-Rejection sampling.
		 *
		 * \copydetails SamplingModule
		 */
		class SamplingModule_PSSACR : public SamplingModule_PDM
		{
		////////////////////////////////
		// Attributes
		protected:
			CompositionRejectionSampler crSampler;

		////////////////////////////////
		// Constructors
		public:
			// Default Constructor
			SamplingModule_PSSACR();
			// Destructor
	virtual	~SamplingModule_PSSACR();

		//////////////////////////////
		// Methods
		protected:
			// Sample reaction index
	virtual bool sample_reaction(pssalib::datamodel::SimulationInfo* ptrSimInfo);
		};
	}
}

#endif /* SAMPLINGMODULE_PSSACR_H_ */
