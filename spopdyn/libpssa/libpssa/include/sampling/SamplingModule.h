/**
 * \file SamplingModule.h
 * \brief Sampling module is responsible for sampling next reaction.
 *
 * \date   13.12.2010
 * \author sanyi
 */

#ifndef SAMPLINGMODULE_H_
#define SAMPLINGMODULE_H_

#include "../typedefs.h"
#include "./CompositionRejectionSampler.h"

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
	namespace sampling
	{
		/**
		 * \class SamplingModule
		 * \brief This class serves as a base class for all sampling schemes used within the library.
		 *
		 * \details Provides the functionality necessary to sample next reaction time and index.
		 * \c SamplingModule class uses the GSL random number generator.
		 * All protected methods of this class may be overloaded, offering a flexible way
		 * to extend its functionality.
		 */
		class SamplingModule
		{
		////////////////////////////////
		// Attributes
		protected:
			//! Pseudo-Random numbers generator from GSL
			gsl_rng * ptrRNG;
			CompositionRejectionSampler crVolumeSampler;

		////////////////////////////////
		// Constructors
		public:
			// Default Constructor
			SamplingModule();

			// Destructor
	virtual ~SamplingModule();

		//////////////////////////////
		// Methods
		protected:
			// Sample next reaction time
	virtual bool sample_time(pssalib::datamodel::SimulationInfo* ptrSimInfo);
	
			// Sample subvolume
	virtual bool sample_volume(pssalib::datamodel::SimulationInfo* ptrSimInfo);

			//! Sample next reaction index
	virtual bool sample_reaction(pssalib::datamodel::SimulationInfo* ptrSimInfo)	=	0;

		public:
			void set_rng_seed(UNSIGNED_INTEGER seed);

			// Get next sample
			bool get_sample(pssalib::datamodel::SimulationInfo* ptrSimInfo);
		};
	}
}

#endif /* GROUPINGMODULE_H_ */
