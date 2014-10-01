/**
 * \file DataModel_DM.h
 * \brief Declares a datatype containing all the datastructures used in the Gillespie's Direct Method.
 *
 * \date   14.12.2010
 * \author sanyi
 */

#ifndef DATAMODEL_DM_H_
#define DATAMODEL_DM_H_

#include "./DataModel.h"

namespace pssalib
{
	namespace datamodel
	{
		////////////////////////////////////////
		// Data model for the Gillespie's direct method
		class DataModel_DM : public DataModel
		{
		public:
			struct Subvolume_DM
			{
				Subvolume_DM();
				~Subvolume_DM();

				//! Array of reaction propensities
				REAL *arPi;
			};

		public:
			DataModel_DM();
			DataModel_DM(DataModel &);
	virtual	~DataModel_DM();
	virtual void Cleanup();

		public:
			Subvolume_DM* arSubvolumes_DM;
		};
	}
}

#endif /* DATAMODEL_DM_H_ */
