/**
 * \file DataModels.cpp
 * \brief Implements \b all classes derived from \c DataModel class.
 *
 * \date   14.12.2010
 * \author sanyi
 */


#include "../../include/PSSA.h"

#include "../../include/datamodel/DataModel.h"
#include "../../include/datamodel/DataModel_DM.h"
#include "../../include/datamodel/DataModel_PDM.h"
#include "../../include/datamodel/DataModel_SPDM.h"
#include "../../include/datamodel/DataModel_PSSACR.h"

namespace pssalib
{
	namespace datamodel
	{
		DataModel::Subvolume::Subvolume()
			: volumeIndex(0)
			, aruN(NULL)
			, dTotalPropensity(0.0)
			, unDiffusionReactions(0)
			, arNeighborVolumes(NULL)
			, ardC(NULL)
			, ardD(NULL)
			, arReactClass(NULL)
			, arDelay(NULL)
			, dH(0.0)
			, dOmega(0.0)
		{
		}

		DataModel::Subvolume::~Subvolume()
		{
			// Clean-up
			if(aruN != NULL)
				delete [] aruN;
			if(arNeighborVolumes != NULL)
				delete [] arNeighborVolumes;
			if(ardC != NULL)
				delete [] ardC;
			if(ardD != NULL)
				delete [] ardD;
			if(arReactClass != NULL)
				delete [] arReactClass;
			if(arDelay != NULL)
				delete [] arDelay;
		}

		////////////////////////////////////////
		// Generic data model class

		DataModel::DataModel()
			: unSpecies(0)
			, unSubvolumes(0)
			, arSubvolumes(NULL)
			, unMiddleVolumeIx(0)
			, dTotalPropensity(0.0)
			, mu(0)
			, sv(0)
			, svDst(0)
		{
		}

		//! Copy constructor
		//! All arrays get detached from source
		//! & are attached to this instance
		DataModel::DataModel(DataModel &d)	:
			strModelName(d.strModelName),
			arU(d.arU),	arUp(d.arUp), arUm(d.arUm),
			unSpecies(d.unSpecies),
			unSubvolumes(d.unSubvolumes),
			arSubvolumes(d.arSubvolumes),
			unMiddleVolumeIx(d.unMiddleVolumeIx)
		{
			d.unSubvolumes = 0;
			d.arSubvolumes = NULL;
			d.strModelName.clear();
		}

		//! Destructor
		DataModel::~DataModel()
		{
			Cleanup();
		}

		void DataModel::Cleanup()
		{
			if(arSubvolumes != NULL)
				delete [] arSubvolumes;
			arSubvolumes = NULL;
		}

		////////////////////////////////////////
		// DM data model class

		DataModel_DM::Subvolume_DM::Subvolume_DM()
			: arPi(NULL)
		{
		}

		DataModel_DM::Subvolume_DM::~Subvolume_DM()
		{
			if(arPi != NULL)
				delete [] arPi;
		}

		//! Default constructor
		DataModel_DM::DataModel_DM()
			: arSubvolumes_DM(NULL)
		{
		}

		//! Copy constructor
		//! All arrays get detached from source
		//! & are attached to this instance
		DataModel_DM::DataModel_DM(DataModel &d)
			: DataModel(d)
			, arSubvolumes_DM(NULL)
		{
		}

		//! Destructor
		DataModel_DM::~DataModel_DM()
		{
			Cleanup();
		}

		void DataModel_DM::Cleanup()
		{
			if(arSubvolumes_DM != NULL)
				delete [] arSubvolumes_DM;
			arSubvolumes_DM = NULL;
			DataModel::Cleanup();
		}

		////////////////////////////////////////
		// PDM data model class
		
		DataModel_PDM::Subvolume_PDM::Subvolume_PDM()
			: arLambda(NULL)
			, arSigma(NULL)
		{
		}
		
		DataModel_PDM::Subvolume_PDM::~Subvolume_PDM()
		{
			if(arLambda != NULL)
				delete [] arLambda;
			if(arSigma != NULL)
				delete [] arSigma;
		}

		//! Default constructor
		DataModel_PDM::DataModel_PDM()
			: arSubvolumes_PDM(NULL)
		{
		}

		//! Copy constructor
		//! All arrays get detached from source
		//! & are attached to this instance
		DataModel_PDM::DataModel_PDM(DataModel &d)
			: DataModel(d)
			, arSubvolumes_PDM(NULL)
		{
		}

		//! Destructor
		DataModel_PDM::~DataModel_PDM()
		{
			Cleanup();
		}

		void DataModel_PDM::Cleanup()
		{
			if(arSubvolumes_PDM != NULL)
				delete [] arSubvolumes_PDM;
			arSubvolumes_PDM = NULL;
			DataModel::Cleanup();
		}

		////////////////////////////////////////
		// SPDM data model class
		
		DataModel_SPDM::Subvolume_SPDM::Subvolume_SPDM()
			: arIdxSigma(NULL)
		{
		}
		
		DataModel_SPDM::Subvolume_SPDM::~Subvolume_SPDM()
		{
			if(arIdxSigma != NULL)
				delete [] arIdxSigma;
		}

		//! Default constructor
		DataModel_SPDM::DataModel_SPDM()
			: arSubvolumes_SPDM(NULL)
			, uI(0)
			, uJ(0)
			, uV(0)
		{
		}

		//! Copy constructor
		//! All arrays get detached from source
		//! & are attached to this instance
		DataModel_SPDM::DataModel_SPDM(DataModel &d)
			: DataModel_PDM(d)
			, arSubvolumes_SPDM(NULL)
			, uI(0)
			, uJ(0)
			, uV(0)
		{
		}

		//! Destructor
		DataModel_SPDM::~DataModel_SPDM()
		{
			Cleanup();
		}

		void DataModel_SPDM::Cleanup()
		{
			if(arSubvolumes_SPDM != NULL)
				delete [] arSubvolumes_SPDM;
			arSubvolumes_SPDM = NULL;
			DataModel_PDM::Cleanup();
		}

		////////////////////////////////////////
		// PSRD-CR data model class

		//! Default constructor
		DataModel_PSSACR::DataModel_PSSACR()
			: arSubvolumes_PSSACR(NULL)
		{
		}

		//! Copy constructor
		//! All arrays get detached from source
		//! & are attached to this instance
		DataModel_PSSACR::DataModel_PSSACR(DataModel &d)
			: DataModel_PDM(d)
			, arSubvolumes_PSSACR(NULL)
		{
		}

		//! Destructor
		DataModel_PSSACR::~DataModel_PSSACR()
		{
			Cleanup();
		}

		void DataModel_PSSACR::Cleanup()
		{
			if(arSubvolumes_PSSACR != NULL)
				delete [] arSubvolumes_PSSACR;
			arSubvolumes_PSSACR = NULL;
			DataModel_PDM::Cleanup();
		}
	}
}
