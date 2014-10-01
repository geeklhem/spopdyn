/**
 * \file  DataModel_PSSACR.h
 * \brief Declares a datatype containing all the datastructures used in the Partial Propensity Direct Method with Composition-Rejection sampling. (Ramaswamy, 2010)
 *
 * \copydetails DataModel.h
 *
 * \date   22.02.2011
 * \author sanyi
 */


#ifndef DATAMODEL_PSSACR_H_
#define DATAMODEL_PSSACR_H_

#include "./DataModel_PDM.h"
#include "./CompositionRejectionSamplerData.h"

namespace pssalib
{
	namespace datamodel
	{
		////////////////////////////////////////
		//! Data model for Partial propensity Direct Method
		//! with composition-rejection sampling
		class DataModel_PSSACR : public DataModel_PDM
		{
		public:
			struct Subvolume_PSSACR
			{
				CompositionRejectionSamplerData crsdSigma;
				std::vector<CompositionRejectionSamplerData> crsdPi;
			};

		public:
			DataModel_PSSACR();
			DataModel_PSSACR(DataModel &);
	virtual ~DataModel_PSSACR();
	virtual void Cleanup();

		public:
			Subvolume_PSSACR* arSubvolumes_PSSACR;
		};
	}
}

#endif /* DATAMODEL_PSSACR_H_ */
