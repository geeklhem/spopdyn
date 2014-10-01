/**
 * \file DataModel_PDM.h
 * \brief Declares a datatype containing all the datastructures used in the Partial Propensity Direct Method. (Ramaswamy, 2009)
 *
 * \date   14.12.2011
 * \author sanyi
 */

#ifndef DATAMODEL_PDM_H_
#define DATAMODEL_PDM_H_

#include "./DataModel.h"

namespace pssalib
{
	namespace datamodel
	{
		////////////////////////////////////////
		// Data model for Partial propensity direct method
		class DataModel_PDM : public DataModel
		{
		public:
			//! Struct for storing the row & column of the dependent propensity
			typedef struct tagPropensityIndex
			{
				UNSIGNED_INTEGER	i, j;

				tagPropensityIndex(UNSIGNED_INTEGER _i,
						UNSIGNED_INTEGER _j) :
							i(_i), j(_j)	{	};

				tagPropensityIndex() :	i(0), j(0)	{	};

				~tagPropensityIndex()	{	};

		inline	bool operator<(const tagPropensityIndex& right) const
				{
					if(i == right.i)
						return (bool)(j < right.j);
					return (bool)(i < right.i);
				};

				// Operator<< for console output
				friend std::ostream & operator<<( std::ostream & output, const tagPropensityIndex & pi)
				{
					output << '(' << pi.i << ',' << pi.j << ')';
					return output;
				};
			} PropensityIndex;

			struct Subvolume_PDM
			{
				Subvolume_PDM();
				~Subvolume_PDM();

				//! Array of arrays of reaction partial propensities.
				MatrixVarLen<REAL>				arPi;

				//! Propensity of each group
				REAL							*arLambda;
				//! Total propensity of each group
				REAL							*arSigma;
			};

		public:
			DataModel_PDM();
			DataModel_PDM(DataModel &);
	virtual ~DataModel_PDM();
	virtual void Cleanup();
			double comb(int n, int coeff);

		public:
			Subvolume_PDM* arSubvolumes_PDM;

			//! Indices of the propensities that need to be updated after
			//! a given reaction has fired.
			MatrixVarLen<PropensityIndex>	arU3;

			//! Look-up table to translate from position in \c arPi to reaction index.
			MatrixVarLen<UNSIGNED_INTEGER>	aruL;
		};
	}
}

#endif /* DATAMODEL_PDM_H_ */
