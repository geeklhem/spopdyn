#pragma once

#include "../typedefs.h"
#include "./PSSACR_Bins.h"

namespace pssalib
{
	namespace datamodel
	{
		class CompositionRejectionSamplerData
		{
		public:
	inline	CompositionRejectionSamplerData()
				: minValue(0.0)
			{}

	inline	CompositionRejectionSamplerData(const CompositionRejectionSamplerData& other)
				: minValue(other.minValue)
				, bins(other.bins)
			{}

	inline	CompositionRejectionSamplerData& operator=(const CompositionRejectionSamplerData& other)
			{
				minValue = other.minValue;
				bins = const_cast<PSSACR_Bins&>(other.bins);
				return *this;
			}

	inline	void updateValue(UNSIGNED_INTEGER bin_no_new, UNSIGNED_INTEGER idx, REAL val)
			{
				bins.updateValue(bin_no_new, idx, val);
			}

		public:
			REAL minValue;
			PSSACR_Bins bins;
		};
	}
}
