/*
 * PSSA.cpp
 *
 *  Created on: 15 груд. 2010
 *      Author: sanyi
 */

#include "../include/PSSA.h"
#include "../include/datamodel/DataModel.h"
#include "../include/datamodel/DataModel_DM.h"
#include "../include/datamodel/DataModel_PDM.h"
#include "../include/datamodel/DataModel_SPDM.h"
#include "../include/datamodel/DataModel_PSSACR.h"

#include "../include/grouping/GroupingModule.h"
#include "../include/grouping/GroupingModule_DM.h"
#include "../include/grouping/GroupingModule_PDM.h"
#include "../include/grouping/GroupingModule_SPDM.h"
#include "../include/grouping/GroupingModule_PSSACR.h"

#include "../include/sampling/SamplingModule.h"
#include "../include/sampling/SamplingModule_DM.h"
#include "../include/sampling/SamplingModule_PDM.h"
#include "../include/sampling/SamplingModule_SPDM.h"
#include "../include/sampling/SamplingModule_PSSACR.h"

#include "../include/update/UpdateModule.h"
#include "../include/update/UpdateModule_DM.h"
#include "../include/update/UpdateModule_PDM.h"
#include "../include/update/UpdateModule_SPDM.h"
#include "../include/update/UpdateModule_PSSACR.h"

#include "../include/datamodel/SimulationInfo.h"

namespace pssalib
{

	PSSA::PSSA()
		: ptrProgrCallback(NULL)
		, ptrErrCallback(NULL)
		, ptrReactionCallback(NULL)
		, ptrProgrCallbackUserData(NULL)
		, ptrErrCallbackUserData(NULL)
		, ptrReactionCallbackUserData(NULL)
		, ptrData(NULL)
		, ptrGrouping(NULL)
		, ptrSampling(NULL)
		, ptrUpdate(NULL)
		, sMethod(M_Invalid)
	{
	}

	PSSA::~PSSA()
	{
		setMethod(M_Invalid);
	}

	void PSSA::SetProgressCallback(FCN_REPORTPROGRESS_CALLBACK fcnProgress, void* user)
	{
		ptrProgrCallback = fcnProgress;
		ptrProgrCallbackUserData = user;
	}

	void PSSA::SetErrorCallback(FCN_REPORTERROR_CALLBACK fcnErr, void* user)
	{
		ptrErrCallback = fcnErr;
		ptrErrCallbackUserData = user;
	}

	void PSSA::SetReactionCallback(FCN_REACTION_CALLBACK fcnReaction, void* user)
	{
		ptrReactionCallback = fcnReaction;
		ptrReactionCallbackUserData = user;
	}
	
	STRING PSSA::GetMethodName(EMethod m)
	{
		switch(m) {
			case M_DM:		return "DM";
			case M_PDM:		return "PDM";
			case M_PSSACR:	return "PSSACR";
			case M_SPDM:	return "SPDM";
			default:		return "Unknown method";
		}
		return "";
	}

	bool PSSA::MakeDir(STRING& path)
	{
#ifdef _WIN32
		if(_mkdir(path.c_str()) != 0)
#else
		if(mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO) != 0)
#endif
		{
			if(errno != EEXIST)
			{
				std::cerr << "mkdir() - could not create output directory!" << std::endl;
				return false;
			}
		}
		return true;
	}

	///////////////////////////////
	// Getter & Setters

	// Sets the method used in simulations
	bool PSSA::setMethod(EMethod sNewMethod)
	{
		if(sMethod != sNewMethod)
		{
			datamodel::DataModel		*oldData = ptrData;
			grouping::GroupingModule	*oldGrouping = ptrGrouping;

			if(ptrSampling != NULL)
				delete ptrSampling;
			ptrSampling = NULL;

			if(ptrUpdate != NULL)
				delete ptrUpdate;
			ptrUpdate = NULL;

			switch(sNewMethod)
			{
			// Unset
			case M_Invalid:
				ptrData = NULL;
				ptrGrouping = NULL;
				break;
			// Gillespie's DM
			case M_DM:
				if(NULL != ptrData)
					ptrData = new datamodel::DataModel_DM(*ptrData);
				else
					ptrData = new datamodel::DataModel_DM();

				if(NULL != ptrGrouping)
					ptrGrouping = new grouping::GroupingModule_DM(*ptrGrouping);
				else
					ptrGrouping = new grouping::GroupingModule_DM();

				ptrSampling	= new sampling::SamplingModule_DM();
				ptrUpdate = new update::UpdateModule_DM();
				break;
			// (Delayed) Partial Propensity Direct Method
			case M_PDM:
				if(NULL != ptrData)
					ptrData = new datamodel::DataModel_PDM(*ptrData);
				else
					ptrData = new datamodel::DataModel_PDM();

				if(NULL != ptrGrouping)
					ptrGrouping = new grouping::GroupingModule_PDM(*ptrGrouping);
				else
					ptrGrouping = new grouping::GroupingModule_PDM();

				ptrSampling	= new sampling::SamplingModule_PDM();
				ptrUpdate = new update::UpdateModule_PDM();
				break;
			// (Delayed) PSSA with Composition-Rejection Sampling
			case M_PSSACR:
				if(NULL != ptrData)
					ptrData = new datamodel::DataModel_PSSACR(*ptrData);
				else
					ptrData = new datamodel::DataModel_PSSACR();

				if(NULL != ptrGrouping)
					ptrGrouping = new grouping::GroupingModule_PSSACR(*ptrGrouping);
				else
					ptrGrouping = new grouping::GroupingModule_PSSACR();

				ptrSampling	= new sampling::SamplingModule_PSSACR();
				ptrUpdate = new update::UpdateModule_PSSACR();
				break;
			// (Delayed) Sorting Partial Propensity Direct Method
			case M_SPDM:
				if(NULL != ptrData)
					ptrData = new datamodel::DataModel_SPDM(*ptrData);
				else
					ptrData = new datamodel::DataModel_SPDM();

				if(NULL != ptrGrouping)
					ptrGrouping = new grouping::GroupingModule_SPDM(*ptrGrouping);
				else
					ptrGrouping = new grouping::GroupingModule_SPDM();

				ptrSampling	= new sampling::SamplingModule_SPDM();
				ptrUpdate = new update::UpdateModule_SPDM();
				break;
			//////////////////////////
			// Illegal parameter value
			default:
				std::cerr << "Error in PSSA::setMethod() : Illegal method identifier '" << sNewMethod << "'.\n";
				return false;
			}

			// Delete old data structures
			if(oldData != NULL)
				delete	oldData;
			if(oldGrouping != NULL)
				delete	oldGrouping;

			// Store new value
			sMethod = sNewMethod;
		}
		return true;
	}

	///////////////////////////////////////////////////
	// Simulation

	bool PSSA::initSimulation(datamodel::SimulationInfo* simInfo)
	{
		//////////////////////////////
		// Check output streams

		if (simInfo->mpi.onlyMaster())
		{
		
			// Status
			bExtStatus = simInfo->ofsStatus.is_open();

			// Debug
			bExtDebug = simInfo->ofsDebug.is_open();
		
			// Initialize the debug file
			if(!bExtDebug)
			{
				std::stringstream ssTemp;
				ssTemp << simInfo->strOutput << "PSSA_debug.txt";
				simInfo->ofsDebug.open(ssTemp.str().c_str(), std::fstream::out);
			}

			// Initialize the status file
			if(!bExtStatus || (NULL == ptrProgrCallback) || (NULL == ptrErrCallback))
			{
				std::stringstream ssTemp;
				ssTemp << simInfo->strOutput << "status.txt";
				simInfo->ofsStatus.open(ssTemp.str().c_str(), std::fstream::out);
			}
		}

		//////////////////////////////
		// Self-checks
		if(!isValid())
		{
			simInfo->reportError("PSSA::initSimulation() - the PSSA has not been initialized properly.");
			return false;
		}

		//////////////////////////////
		// Check the parameters

		// Attach SimulationInfo object to this instance
		simInfo->attachPSSA(this);

		// Check if model was loaded
		if(!simInfo->isDataLoaded())
		{
			simInfo->reportError("PSSA::initSimulation() - no model was loaded.");
			return false;
		}

		// Check input params
		if(!simInfo->isValid())
		{
			simInfo->reportError("PSSA::initSimulation() - simulation parameters are invalid.");
			return false;
		}

		// Reset timers
		simInfo->dTimeCheckpoint = 0.0;
		simInfo->dTimeSimulation = 0.0;

		return true;

	}

	bool PSSA::deinitSimulation(datamodel::SimulationInfo*	simInfo)
	{
	  
		if (simInfo->mpi.onlyMaster())
		{
			// Close output streams
			if(!bExtStatus)
				simInfo->ofsStatus.close();

			// Close the debug output stream
			if(!bExtDebug)
				simInfo->ofsDebug.close();
		}
			
		// Dettach SimulationInfo object from this instance
		simInfo->detachPSSA();

		for (std::vector<INTEGER*>::iterator i = initial_population.begin(); i != initial_population.end(); ++i)
			delete [] *i;
		initial_population.clear();

		return true;
	}

	bool PSSA::setupForSampling(datamodel::SimulationInfo* simInfo, bool initSpeciesArray)
	{
		if(!initSimulation(simInfo))
			return false;

		//////////////////////////////
		// Initialize the data structures
		if(!ptrGrouping->preinitialize(simInfo))
		{
			// Failed, report & exit
			std::stringstream ssTemp;
			ssTemp << "Error : Failed to initialize data structures." << std::endl;
			simInfo->reportError(ssTemp.str());
			return false;
		}
		
		if (initSpeciesArray)
			initializeSpeciesArray(simInfo);

		// Store the initial population
		for (std::vector<INTEGER*>::iterator i = initial_population.begin(); i != initial_population.end(); ++i)
			delete [] *i;
		initial_population.clear();
		initial_population.resize(ptrData->unSubvolumes);
		for (UNSIGNED_INTEGER i = 0; i < ptrData->unSubvolumes; i++) {
			initial_population[i] = new INTEGER[ptrData->unSpecies];
			memcpy(initial_population[i], ptrData->arSubvolumes[i].aruN, sizeof(INTEGER)*ptrData->unSpecies);
		}

		return true;
	}

	bool PSSA::doSampleLoop(datamodel::SimulationInfo* simInfo, UNSIGNED_INTEGER unSamples, bool OutputTrajectory, std::ofstream& ofsTiming, std::ofstream& ofsFinalVals,
		std::vector< UNSIGNED_INTEGER >* final_data_points)
	{ 
		UNSIGNED_INTEGER n = 0;
		UNSIGNED_INTEGER n_it = 0;
		simInfo->mpi.init(unSamples);
		while (simInfo->mpi.spread(n))
		{
//		for(UNSIGNED_INTEGER n = 0; n < unSamples; n++)
//		{
			if (simInfo->bVerbose)
				std::cout << "Sample " << n  << std::endl;

			if (OutputTrajectory && 0.0 == simInfo->dTimeStart)
			{
				// Output the initial population to stream
				for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; si++) {
					if (si != 0) simInfo->ssOutput << "\t";
					for(UNSIGNED_INTEGER i = 1; i < ptrData->unSpecies; i++) {
						if (i != 1) simInfo->ssOutput << ",";
						simInfo->ssOutput << ptrData->arSubvolumes[si].aruN[i];
					}
				}
				simInfo->ssOutput << std::endl;
				
				// Vtk output
				
				if (simInfo->vtkOutName.size() != 0)
				{
					std::stringstream rank_str;
					rank_str << "_" << simInfo->mpi.getRank() << "_";
					simInfo->vtkProcess.WriteVtkFile(simInfo->vtkOutName + rank_str.str() ,simInfo,0);
				}
			}

			// Initialize data structures
			if (!ptrGrouping->initialize(ptrData))
			{
				// Failed, report & exit
				std::stringstream ssTemp;
				ssTemp << "Error : Failed to initialize data structures." << std::endl;
				simInfo->reportError(ssTemp.str());
				return false;
			}
			ptrGrouping->postinitialize(ptrData);

			// Zero times
			REAL tTrial = 0.0;
			UNSIGNED_INTEGER unReactions = 0;

			// Report progress of the simulation
			simInfo->reportProgress(n, unSamples);

			// Start timing
			simInfo->beginTrial();

			/////////////////////////////////
			// Run the internal loop
			while(simInfo->isRunning())
			{ 
				bool bResult = ptrSampling->get_sample(simInfo);

				if (OutputTrajectory)
					simInfo->doOutput();

				if(bResult)
				{
					bResult = ptrUpdate->do_update(simInfo);
					++unReactions;
				}
				
				if(bResult)
				{
					if (ptrReactionCallback) {
						ptrReactionCallback(ptrData, simInfo->dTimeSimulation, ptrReactionCallbackUserData);
					}
				}
				else
				{
					if(!simInfo->isRunning())
						break;
					else
						return false;
				}
				
				int totsp = 0;
				
				for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; si++)
					totsp += ptrData->arSubvolumes[si].aruN[2];
			}

			// End timing
			tTrial = simInfo->endTrial();

			// Report progress of the simulation to a file
			// TODO Merge this with endTrial();
			simInfo->reportProgress(n + 1, unSamples);
			
			// Output the timing information
			if (ofsTiming.is_open())
			{
				ofsTiming << std::setw(14) << tTrial << std::setw(14) << unReactions << std::setw(14) << tTrial / (REAL) unReactions << std::endl;
			}

			// Output the population at final time point
			if (ofsFinalVals.is_open())
			{
				for(UNSIGNED_INTEGER si = 0; si < ptrData->unSubvolumes; si++) {
					if (si != 0) simInfo->ssOutput << "\t";
					for(UNSIGNED_INTEGER i = 1; i < ptrData->unSpecies; i++) {
						if (i != 1) simInfo->ssOutput << ",";
						simInfo->ssOutput << ptrData->arSubvolumes[si].aruN[i];
					}
				}
			}

			// Store the population
			
			if (final_data_points && datamodel::SimulationInfo::ST_NONE != simInfo->stStatistics)
			{
				for(UNSIGNED_INTEGER i = 0; i < arSpecies.size(); i++)
				{
					int idx_sp = arSpecies[i];
					(*final_data_points)[i+n_it*arSpecies.size()] = ptrData->arSubvolumes[0].aruN[idx_sp];
					for(UNSIGNED_INTEGER si = 1; si < ptrData->unSubvolumes; si++)
					{
						(*final_data_points)[n_it*arSpecies.size()+i] += ptrData->arSubvolumes[si].aruN[idx_sp];
					}
				}
			}
			
			if (OutputTrajectory)
			{
				// Output the final result
				std::stringstream ssTemp;
				ssTemp << simInfo->strOutput << "PSSA_trajectory_" << n << "_" << simInfo->mpi.getRank() << ".txt";

				std::ofstream ofsOutput;
				ofsOutput.open(ssTemp.str().c_str(), std::fstream::out);
				ofsOutput << simInfo->ssOutput.rdbuf();
				ofsOutput.close();

				// Clear streams
				simInfo->ssOutput.str("");
			}

			// Reset population
			for (UNSIGNED_INTEGER i = 0; i < ptrData->unSubvolumes; i++) {
				memcpy(ptrData->arSubvolumes[i].aruN, initial_population[i], sizeof(INTEGER)*ptrData->unSpecies);
			}
			
			n_it++;
		}
		
		if (final_data_points != NULL)
			simInfo->mpi.collect_result(&(*final_data_points)[0],&(*final_data_points)[0],unSamples,arSpecies.size()*sizeof(UNSIGNED_INTEGER));
		
		return true;
	}

	void PSSA::initializeSpeciesArray(datamodel::SimulationInfo* simInfo)
	{
		if((simInfo->arSpeciesIds.size() == 1)&&(simInfo->arSpeciesIds[0].compare("all species") == 0))
		{
			for(UNSIGNED_INTEGER i = 1; i < ptrData->unSpecies; ++i)
				arSpecies.push_back(i);
		}
		else
		{
			grouping::GroupingModule::MAP_SP2IDX::iterator	it;
			for(UNSIGNED_INTEGER i = 0; i < simInfo->arSpeciesIds.size(); ++i)
			{
				if (simInfo->mpi.onlyMaster())
					simInfo->ofsDebug << "Element #" << i << " : ";

				it = ptrGrouping->mapSp2Idx.find(simInfo->arSpeciesIds[i]);
				if(it != ptrGrouping->mapSp2Idx.end())
				{
					if (simInfo->mpi.onlyMaster())
						simInfo->ofsDebug << "adding species with identifier '" << simInfo->arSpeciesIds[i] << "' as species #" << it->second << std::endl;
					arSpecies.push_back(it->second);
				}
				else
				{
					std::stringstream ssTemp;
					if (simInfo->mpi.onlyMaster())
						ssTemp << "PSSA::run_hist() - species identifier '" << simInfo->arSpeciesIds[i] << "' not found in model; position " << i;
					simInfo->reportError(ssTemp.str());
				}
			}
		}
	}

	/**
   * This function samples given number of trajectories (\c simInfo->arSamples[0]) and outputs them to a series of files, 
   * starting at \a time \c = \c simInfo->dTimeStart to \a time \c = \c simInfo->dTimeEnd seconds and saving 
   * the output to \c simInfo->strOutput every \c simInfo->dTimeStep seconds. If \c simInfo->stStatistics is set 
   * to \a ST_SINGLE, statistical information, such as mean trajectory and standard deviation from mean trajectory are output to a special file.
   *
   * @param simInfo datamodel::SimulationInfo* Simulation information object associated with this run.
   *
   * @return logical \a "true" if the all trials finish successfully,
   * \a "false" otherwise.
   *
   */
	bool PSSA::run_avg(datamodel::SimulationInfo*	simInfo)
	{
		std::ofstream ofsOutput, ofsTiming, ofsFinalVals;

		simInfo->mpi.redirectOut();
		
		if (!setupForSampling(simInfo, false))
			return false;
		
		UNSIGNED_INTEGER unSamples = simInfo->arSamples[0];
		if (!doSampleLoop(simInfo, unSamples, true, ofsTiming, ofsFinalVals, NULL))
		{
			deinitSimulation(simInfo);
			return false;
		}
		
		// Species names
		{
			std::stringstream ssTemp;
			ssTemp << simInfo->strOutput << "PSSA_names.txt";

			ofsOutput.open(ssTemp.str().c_str(), std::fstream::out);
			for(UNSIGNED_INTEGER i = 0; i < ptrGrouping->arIdx2Sp.size(); ++i)
				ofsOutput << ptrGrouping->arIdx2Sp[i] << std::endl;
			ofsOutput.close();
		}

#if BT_TODO // This will need some work; not sure how useful it is in the reaction-diffusion case.
		// Statistics
		if(datamodel::SimulationInfo::ST_SINGLE == simInfo->stStatistics)
		{
			// Temporary arrays storing mean values & standard deviation.,
			// timing vars & others
			REAL				*mean, *stddev;
			UNSIGNED_INTEGER	tempi = 0, N = 0;
			// Temporary streams
			std::stringstream	ssMean, ssStdDev;
			// Input from file
			std::ifstream		ifsInput;

			// Allocate memory
			mean	= new REAL[ptrData->unSpecies - 1];
			stddev	= new REAL[ptrData->unSpecies - 1];

			// Reset timing
			simInfo->dTimeCheckpoint = simInfo->dTimeStart;

			// TODO Code for generating statistical info
			while(simInfo->dTimeCheckpoint < simInfo->dTimeEnd)
			{
				memset(mean, 0, sizeof(REAL)*(ptrData->unSpecies - 1));
				memset(stddev, 0, sizeof(REAL)*(ptrData->unSpecies - 1));

				// Compute mean
				for(UNSIGNED_INTEGER n = 0; n < unSamples; n++)
				{
					// Read the population at this point of time in trial n
					std::stringstream ssTemp;
					ssTemp << simInfo->strOutput << "PSSA_trajectory_" << n << ".txt";
					ifsInput.open(ssTemp.str().c_str(), std::fstream::out);

					// Skip other rows
					for(UNSIGNED_INTEGER j = 0; j < N; j++)
					{
						STRING tmp;
						getline(ifsInput, tmp);
					}

					for(UNSIGNED_INTEGER j = 0; j < (ptrData->unSpecies - 1); j++)
					{
						ifsInput >> tempi;
						mean[j] += tempi;
					}
					ifsInput.close();
				}

				// Output mean
				for(UNSIGNED_INTEGER j = 0; j < (ptrData->unSpecies - 1); j++)
				{
					mean[j]	= mean[j] / ((REAL) unSamples);
					ssMean << std::setw(12) << mean[j];
				}
				ssMean << std::endl;

				// Compute stddev
				for(UNSIGNED_INTEGER n = 0; n < unSamples; n++)
				{
					// Read the population at this point of time in trial n
					std::stringstream ssTemp;
					ssTemp << simInfo->strOutput << "PSSA_trajectory_" << n << ".txt";
					ifsInput.open(ssTemp.str().c_str(), std::fstream::out);

					// Skip other rows
					for(UNSIGNED_INTEGER j = 0; j < N; j++)
					{
						STRING tmp;
						getline(ifsInput, tmp);
					}

					for(UNSIGNED_INTEGER j = 0; j < (ptrData->unSpecies - 1); j++)
					{
						ifsInput >> tempi;
						stddev[j] += pow(((REAL) tempi) - mean[j], 2);
					}
					ifsInput.close();
				}

				// Output stddev
				for(UNSIGNED_INTEGER j = 0; j < (ptrData->unSpecies - 1); j++)
				{
					stddev[j]	= sqrt(stddev[j] / ((REAL) unSamples));
					ssStdDev << std::setw(12) << stddev[j];
				}
				ssStdDev << std::endl;

				// Advance the indexes
				++N;	simInfo->dTimeCheckpoint += simInfo->dTimeStep;
			}

			// clean-up
			delete	[]	mean;
			delete	[]	stddev;

			// Output the mean trajectory
			{
				std::stringstream ssTemp;
				ssTemp << simInfo->strOutput << "PSSA_mean_trajectory.txt";

				ofsOutput.open(ssTemp.str().c_str(), std::fstream::out);
				ofsOutput << ssMean.rdbuf();
				ofsOutput.close();
			}

			// Output the mean trajectory
			{
				std::stringstream ssTemp;
				ssTemp << simInfo->strOutput << "PSSA_std_dev.txt";

				ofsOutput.open(ssTemp.str().c_str(), std::fstream::out);
				ofsOutput << ssStdDev.rdbuf();
				ofsOutput.close();
			}
		}
#endif

		//
		// Clean-up
		deinitSimulation(simInfo);
		simInfo->mpi.closeOut();
		
		return true;
	}

	/**
   * This function simulates given in \c simInfo->arSamples number of trajectories and outputs the population at \a time \c = \c simInfo->dTimeEnd 
   * to a file named \c "PSSA_final_values.txt". For each element in \c arTrials a separate output is produced and the data from previous trials (if any) is
   * reused. At the same time it computes the amount of time spent one each run and the number of reactions fired and stores them to \c "PSSA_probability_timing.txt". 
   * Depending on the value of \c simInfo->stStatistics, it can produce the one-dimensional histograms for each of the species (\a ST_SINGLE) or cumulative histogram 
   * of all of the species (\a ST_MULTI). Only species in \c simInfo->arSpeciesIds are considered.
   *
   * @param simInfo datamodel::SimulationInfo* Simulation information object associated with this run.
   *
   * @return logical \a "true" if the all trials finish successfully,
   * \a "false" otherwise.
   *
   */
	bool PSSA::run_hist(pssalib::datamodel::SimulationInfo *simInfo)
	{
		std::ofstream ofsTiming, ofsFinalVals, ofsOutput;
		
		simInfo->mpi.redirectOut();
		
		if (!setupForSampling(simInfo, true))
			return false;

		// Allocate memory.
		UNSIGNED_INTEGER unSamplesMax = (*(simInfo->arSamples.rbegin()));
		std::vector<UNSIGNED_INTEGER> final_data_points;
		if((datamodel::SimulationInfo::ST_NONE != simInfo->stStatistics) && (arSpecies.size() > 0))
		{
			if (simInfo->mpi.onlyMaster())
				final_data_points.resize(unSamplesMax*arSpecies.size());
			else
				final_data_points.resize(simInfo->mpi.dataChuckSize(unSamplesMax*arSpecies.size()));
		}

		//////////////////////////////
		// Outer samples loop
		for(UNSIGNED_INTEGER i_s = 0, n = 0; i_s < simInfo->arSamples.size(); i_s++)
		{
			// Get current number of samples
			UNSIGNED_INTEGER unSamples = simInfo->arSamples[i_s];

			// Setup output path
			STRING strOutputPath;
			{
				std::stringstream ssTemp;
				ssTemp << simInfo->strOutput << unSamples << '/';
				strOutputPath = ssTemp.str();
			}
			MakeDir(strOutputPath);

			//////////////////////////////
			// Setup the output

			// Timing
			{
				std::stringstream ssTemp;
				ssTemp << strOutputPath << "PSSA_probability_timing.txt";
				ofsTiming.open(ssTemp.str().c_str(), std::fstream::out);
			}

			// Final values
			{
				std::stringstream ssTemp;
				ssTemp << strOutputPath << "PSSA_final_values.txt";
				ofsFinalVals.open(ssTemp.str().c_str(), std::fstream::out);
			}

			if (!doSampleLoop(simInfo, unSamples, false, ofsTiming, ofsFinalVals, &final_data_points))
			{
				deinitSimulation(simInfo);
				return false;
			}

			// Close output streams
			ofsTiming.close();
			ofsFinalVals.close();

			// Clear streams
			simInfo->ssOutput.str("");
			if(datamodel::SimulationInfo::SimulationInfo::ST_NONE != simInfo->stStatistics && simInfo->mpi.onlyMaster() )
			{
				// Output variables
				REAL										unSamplesInv =	1.0 / (REAL) unSamples;
				std::vector< std::vector<UNSIGNED_INTEGER> >::iterator 	it, curr_it;
				UNSIGNED_INTEGER							unOccur;

				//////////////////////////////
				// Output loop
				
				std::vector< std::vector<UNSIGNED_INTEGER> > final_data_points_sort;

				// Fill the data for sort process
				
				final_data_points_sort.resize(unSamples);

				for (unsigned int j = 0 ; j < unSamples ; j++)
				{
					for (unsigned int i = 0 ; i < arSpecies.size() ; i++)
					{
						final_data_points_sort[j].push_back(final_data_points[j*arSpecies.size()+i]);
					}
				}
				
				for(UNSIGNED_INTEGER i_o = 0; i_o < arSpecies.size(); i_o++)
				{
					unOccur = ((n / 2) > 0 ? (n / 2) : 1);
					
					std::vector< REAL > P;
					std::vector< std::vector<UNSIGNED_INTEGER> > unique_data_points;
					
					P.reserve(unOccur);
					unique_data_points.reserve(unOccur);

					// Set DataSorter mode
					DataSorter sorter;
					if(datamodel::SimulationInfo::ST_SINGLE == simInfo->stStatistics)
						sorter.SetIndex(i_o);
					else
						sorter.UnsetIndex();

					// Sort the result
					std::sort(final_data_points_sort.begin(), final_data_points_sort.end(), sorter);

					// Calculate the Cumulative Probability Density function
					it = final_data_points_sort.begin();	curr_it = it;	unOccur = 0;
					while(final_data_points_sort.end() != it)
					{
						if(!sorter.equals((*curr_it),(*it)))
						{
							unique_data_points.push_back(*curr_it);
							P.push_back( (REAL)unOccur * unSamplesInv );
							unOccur = 0;
							curr_it = it;
						}
						++unOccur;
						++it;
					}

					// Store the last element
					
					unique_data_points.push_back(*curr_it);
					P.push_back( (REAL)unOccur * unSamplesInv );

					// Dump probability distribution to a file
					unOccur = unique_data_points.size();
					for(UNSIGNED_INTEGER j = 0; j < unOccur; j++)
						simInfo->ssOutput << sorter.getString(unique_data_points[j]) << std::setw(1) << ' ' << std::setw(12) << std::setprecision(9) << P[j] << std::endl;
					simInfo->ofsDebug << "Output for species # " << arSpecies[i_o] << std::flush << " with id '" << ptrGrouping->arIdx2Sp[arSpecies[i_o] - 1] << "'" << std::endl << std::flush;

					std::stringstream ssTemp;
					{
						ssTemp << strOutputPath << "PSSA_probability";
						if(datamodel::SimulationInfo::ST_MULTI != simInfo->stStatistics)
							ssTemp << "_" << ptrGrouping->arIdx2Sp[arSpecies[i_o] - 1];
						ssTemp << ".txt";

						ofsOutput.open(ssTemp.str().c_str(), std::fstream::out);
						ofsOutput << simInfo->ssOutput.rdbuf();
						ofsOutput.close();
					}

					// Clear streams
					simInfo->ssOutput.str("");

					// Clear arrays
					unique_data_points.clear();
					P.clear();

					// Stop if it should be multidimensional
					if(datamodel::SimulationInfo::ST_MULTI == simInfo->stStatistics)
						break;
					
				
					if (simInfo->plotVisualize == true)
						simInfo->plot.GraphHist(ssTemp.str());
				}
				// End of output loop
				//////////////////////////////
			}
		}
		// Outer samples loop
		//////////////////////////////

		//
		// Clean-up
		deinitSimulation(simInfo);
		simInfo->mpi.closeOut();
		
		return true;
	}

	///////////////////////////////////////////////////////////////
	// Global functions specific for PSSAlib

	REAL	getHmu(UNSIGNED_INTEGER	n, UNSIGNED_INTEGER	m)
	{
		if(m > 1)
			return	((REAL)n * getHmu(n - 1, m - 1)) / ((REAL) m);
		else
			return	n;
	}


	UNSIGNED_INTEGER factorial(UNSIGNED_INTEGER n)
	{
	   UNSIGNED_INTEGER	result = 1;

	   for ( ; n > 1; n-- )
	      result *= n;

	   return result;
	}
}
