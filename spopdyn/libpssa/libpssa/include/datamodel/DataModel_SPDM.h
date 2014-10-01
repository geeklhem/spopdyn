/**
 * \file DataModel_SPDM.h
 * \brief Declares a datatype containing all the datastructures used in the Sorting Partial Propensity Direct Method. (Ramaswamy, 2009)
 *
 * \date   14.05.2011
 * \author sanyi
 */


#ifndef DATAMODEL_SPDM_H_
#define DATAMODEL_SPDM_H_

#include "./DataModel_PDM.h"

namespace pssalib
{
	namespace datamodel
	{
		////////////////////////////////////////
		// Data model for Partial propensity direct method
		class DataModel_SPDM : public DataModel_PDM
		{
		public:
			struct Subvolume_SPDM
			{
				Subvolume_SPDM();
				~Subvolume_SPDM();

				//! Indexes for sampling using PDM structures
				UNSIGNED_INTEGER				*arIdxSigma;
				//! Indexes for sampling using PDM structures
				MatrixVarLen<UNSIGNED_INTEGER>	arIdxPi;
			};

		public:
			DataModel_SPDM();
			DataModel_SPDM(DataModel &);
	virtual	~DataModel_SPDM();
	virtual void Cleanup();

		public:
			Subvolume_SPDM* arSubvolumes_SPDM;

			// Position
			UNSIGNED_INTEGER uI, uJ, uV;
		};
	}
}

#endif /* DATAMODEL_SPDM_H_ */
