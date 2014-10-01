#pragma once

#include "../typedefs.h"
#include "../datamodel/CompositionRejectionSamplerData.h"

namespace pssalib
{
	namespace sampling
	{
		class CompositionRejectionSampler
		{
		public:
			bool Sample(const pssalib::datamodel::CompositionRejectionSamplerData* ptrData, gsl_rng* ptrRNG, const REAL scale, UNSIGNED_INTEGER& outI, REAL& outR);
		};
	}
}
