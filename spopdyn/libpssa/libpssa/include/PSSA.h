/**
 * \file   PSSA.h
 * \brief  PSSA library main include header
 *
 * This file contains the definition of the API.
 * It provides access to main functions of the library:
 * - creating a simulation engine;
 * - simulating using the "run"-mode: outputs trajectories
 * - simulating using the "stat"-mode: outputs PDFs.
 *
 * \date   15.12.2010
 * \author sanyi
 */

#ifndef PSSA_H_
#define PSSA_H_

#include "typedefs.h"

// Forward declarations
namespace pssalib
{
	namespace datamodel
	{
		class DataModel;
		class SimulationInfo;

		// Class Factories
		class SimulationInfoClassFactory;
	}

	namespace grouping
	{
		class GroupingModule;
	}

	namespace sampling
	{
		class SamplingModule;
	}

	namespace update
	{
		class UpdateModule;
	}
}

namespace pssalib
{
	///////////////////////////////////
	//! Partial propensity stochastic reaction-diffusion simulation library
	/*!
	 * Main simulation engine class.
	 */
	class PSSA
	{
	public:
		enum EMethod
		{
			//! Invalid Method
			M_Invalid = -1,
			//! Direct Method
			M_DM = 0,
			//! Partial Propensity Direct Method
			M_PDM = 1,
			//! Partial Propensity Sample Rejection
			M_PSSACR = 2,
			//! Sorting Partial Propensity with Sample Rejection
			M_SPDM = 3
		};

		PSSA();
		~PSSA();

		void SetProgressCallback(FCN_REPORTPROGRESS_CALLBACK fcnProgress, void* user);
		void SetErrorCallback(FCN_REPORTERROR_CALLBACK fcnErr, void* user);
		void SetReactionCallback(FCN_REACTION_CALLBACK fcnReaction, void* user);
		
		static STRING GetMethodName(EMethod m);
		static bool MakeDir(STRING& path);
		
		//! Sample given number of trajectories, output them and compute statistics.
		bool run_avg(datamodel::SimulationInfo*	simInfo);
		//! Compute histogram of species at a given time point over a given number of trials.
		bool run_hist(datamodel::SimulationInfo*	simInfo);

	////////////////////////////////////////////////////////////////////////////
	// Class members
	private:
		FCN_REPORTPROGRESS_CALLBACK				ptrProgrCallback;
		FCN_REPORTERROR_CALLBACK				ptrErrCallback;
		FCN_REACTION_CALLBACK					ptrReactionCallback;
		void*									ptrProgrCallbackUserData;
		void*									ptrErrCallbackUserData;
		void*									ptrReactionCallbackUserData;

		datamodel::DataModel*					ptrData;
		grouping::GroupingModule*				ptrGrouping;
		sampling::SamplingModule*				ptrSampling;
		update::UpdateModule*					ptrUpdate;

		//! ID of the method used
		EMethod									sMethod;

		//! String representing last error
		STRING									strError;

		// Temporary variables
		std::vector<UNSIGNED_INTEGER> 			arSpecies;
		std::vector<INTEGER*>					initial_population;

		// Friend classes definition
		friend class datamodel::SimulationInfo;

	protected:
		// Temporary variables (output streams)
		bool				bExtStatus;
		bool				bExtDebug;

	////////////////////////////////////////////////////////////////////////////
	// Internal helpers
	protected:

		//! Run self-checks & prepare the simInfo object for simulations
		bool initSimulation(datamodel::SimulationInfo* simInfo);

		//! Free all resources
		bool deinitSimulation(datamodel::SimulationInfo* simInfo);

		bool setupForSampling(datamodel::SimulationInfo* simInfo, bool initSpeciesArray);
		bool doSampleLoop(datamodel::SimulationInfo* simInfo, UNSIGNED_INTEGER unSamples, bool OutputTrajectory, std::ofstream& ofsTiming, std::ofstream& ofsFinalVals,
			std::vector<UNSIGNED_INTEGER > * final_data_points);
		void initializeSpeciesArray(datamodel::SimulationInfo* simInfo);

	////////////////////////////////////////////////////////////////////////////
	// Getters & Setters
	public:
		//! \internal Set the value of error string
		void						reportError(STRING sError)
		{
			strError = sError;
		};

		//! Return a string containing the description of last error occured
		STRING						getLastError()	const
		{
			return strError;
		};

		//! Returns the method used for simulations
		inline	EMethod				getMethod() const	{	return sMethod;	};

		//! Sets the method used in simulations
		bool						setMethod(EMethod sNewMethod);

		/*
		 * Getters for modules
		 */

		//! Returns a pointer to a DataModel object asociated with current instance of the engine
		inline	STRING				getModelName()			const;

		//! Returns number of species in the model
		inline	UNSIGNED_INTEGER	getNumSpecies()			const;

		//! Returns number of reactions in the model
		inline	UNSIGNED_INTEGER	getNumReactions()		const;

		//! Returns number of queued reactions in the model
		inline	UNSIGNED_INTEGER	getNumQueuedReactions()	const;

		//! Returns true if the engine was properly initialized
		inline	bool				isValid()	const
		{
			if(sMethod != (SHORT)(-1))
				return (ptrData != NULL)&&(ptrGrouping != NULL)&&(ptrSampling != NULL)&&(ptrUpdate != NULL);
			return false;
		}
	};
}

// For the inline method (see below)
#include "./datamodel/DataModel.h"

namespace pssalib
{
	inline	STRING				PSSA::getModelName()			const
	{
		return	ptrData->strModelName;
	}

	inline	UNSIGNED_INTEGER	PSSA::getNumSpecies()			const
	{
		return	ptrData->unSpecies;
	}
}


// To let the user use the SimulationInfo objects
// just by including one file
#include "./datamodel/SimulationInfo.h"

#endif
