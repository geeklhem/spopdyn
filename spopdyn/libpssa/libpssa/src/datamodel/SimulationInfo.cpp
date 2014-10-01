#include "../../include/PSSA.h"
#include "../../include/datamodel/DataModel.h"
#include "../../include/update/UpdateModule.h"
#include "../../include/datamodel/SimulationInfo.h"

LIBSBML_CPP_NAMESPACE_USE

namespace pssalib
{
	namespace datamodel
	{
		/////////////////////////////////
		// Constructors

		//! Constructor
		SimulationInfo::SimulationInfo()
			: ptrPSSA(NULL)
			, bOverrideUpdate(false)
			, dTimeCheckpoint(0.0)
			, dTimeStep(0.0)
			, dTimeEnd(0.0)
			, dTimeSimulation(0.0)
			, stStatistics(ST_NONE)
			, eVolumeType(VT_Invalid)
			, eBoundaryConditions(BC_Invalid)
			, unNumGridPoints(1)
			, dOmega(0.0)
			, eInitialPop(IP_Invalid)
			, uUpdateMode(0)
			, bDisableOutput(false)
			, bVerbose(false)
			, eInitialUserDefined(NULL)
			, frame(0)
			, vtkProcess()
			, vtkOutName()
		{
#ifdef __MACH__
			mach_timebase_info_data_t    sTimebaseInfo;
			mach_timebase_info(&sTimebaseInfo);
			tickRatio = ((REAL_EXT)  sTimebaseInfo.numer) / ((REAL_EXT) sTimebaseInfo.denom) / ((REAL_EXT) 1e9);
#elif defined(_WIN32)
			LARGE_INTEGER frequency;
			QueryPerformanceFrequency(&frequency);
			tickRatio = (REAL_EXT)frequency.QuadPart;
#endif
		}

		//! Copy constructor
		SimulationInfo::SimulationInfo(SimulationInfo & right)
			: ptrPSSA(right.ptrPSSA)
			, bOverrideUpdate(right.bOverrideUpdate)
			, dTimeCheckpoint(right.dTimeCheckpoint)
			, dTimeStep(right.dTimeEnd)
			, dTimeEnd(right.dTimeEnd)
			, dTimeSimulation(right.dTimeSimulation)
			, stStatistics(right.stStatistics)
			, eVolumeType(right.eVolumeType)
			, eBoundaryConditions(right.eBoundaryConditions)
			, unNumGridPoints(right.unNumGridPoints)
			, dOmega(right.dOmega)
			, eInitialPop(right.eInitialPop)
			, uUpdateMode(right.uUpdateMode)
			, bDisableOutput(right.bDisableOutput)
			, bVerbose(right.bVerbose)
			, eInitialUserDefined(right.eInitialUserDefined)
			, frame(0)
			, vtkProcess()
			, vtkOutName()
		{
			ssOutput.str("");
			ssOutput << right.ssOutput.str();
#ifdef __MACH__
			mach_timebase_info_data_t    sTimebaseInfo;
			mach_timebase_info(&sTimebaseInfo);
			tickRatio = ((REAL_EXT)  sTimebaseInfo.numer) / ((REAL_EXT) sTimebaseInfo.denom) / ((REAL_EXT) 1e9);
#elif defined(_WIN32)
			LARGE_INTEGER frequency;
			QueryPerformanceFrequency(&frequency);
			tickRatio = (REAL_EXT)frequency.QuadPart;
#endif
		}

		//! Destructor
		SimulationInfo::~SimulationInfo()
		{
		}

		//////////////////////////////
		// Methods

		// Check this instance
		bool SimulationInfo::isValid()	
		{
			// Check if this instance is attached
			if(NULL == ptrPSSA)
			{
				reportError("SimulationInfo::isValid() - Error: method should not be called on a detached instance.");
				return false;
			}

			// Check output path
			if(!chkPath(strOutput))
			{
				reportError("SimulationInfo::isValid() - Error: invalid output path.");
				return false;
			}

			//
			// Check if number of trial is valid

			// Sort trials array
			std::sort(arSamples.begin(),arSamples.end());

			// Check if the biggest one is greater than zero
			if(0 == (*(arSamples.rbegin())))
			{
				reportError("SimulationInfo::isValid() - Error: number of trials must be a positive integer( > 0).");
				return false;
			}

			if (dOmega < 0) {
				std::cerr << "SimulationInfo::isValid: Negative omega: " << dOmega << std::endl;
				return false;
			}
			if (unNumGridPoints == 0) {
				std::cerr << "SimulationInfo::isValid: Zero grid points" << std::endl;
				return false;
			}
			if (eVolumeType == VT_Invalid) {
				std::cerr << "SimulationInfo::isValid: Invalid volume type" << std::endl;
				return false;
			}
			if (eBoundaryConditions == BC_Invalid) {
				std::cerr << "SimulationInfo::isValid: Invalid boundary conditions" << std::endl;
				return false;
			}
			if (eInitialPop == IP_Invalid) {
				std::cerr << "SimulationInfo::isValid: Invalid initial population type" << std::endl;
				return false;
			}

			return true;
		}

		// load a new model file
		bool SimulationInfo::readSBMLFile( STRING strInput )
		{
			// Try to load SBML file
			ptrSBMLDocument	= readSBML(strInput.c_str());

			if (ptrSBMLDocument->getNumErrors() > 0)
			{
				ptrSBMLDocument->printErrors(std::cerr);
				return false;
			}
			return true;
		}

		//! Initializes timing variables
		void SimulationInfo::beginTrial()
		{
			// initialize the offset for output
			if(dTimeStart > 0.0)
				dTimeCheckpoint = dTimeStart - dTimeStep;
#ifdef __linux__
			clock_gettime(CLOCK_MONOTONIC, &trialStart);
#elif defined(__MACH__)
			trialStart = mach_absolute_time();
#elif defined(_WIN32)
			QueryPerformanceCounter(&trialStart);
#endif
		}

		//! Measures time spent on last trial (in seconds) and resets the timers.
		REAL_EXT SimulationInfo::endTrial()
		{
#ifdef __linux__
			// get hi res time
			clock_gettime(CLOCK_MONOTONIC, &trialEnd);
#elif defined(__MACH__)
			// Get absolute time
			trialEnd = mach_absolute_time();
#elif defined(_WIN32)
			QueryPerformanceCounter(&trialEnd);
#endif
			// Output the last change (if any)
			doOutput();

			// reset the timing
			dTimeSimulation = 0.0;
			dTimeCheckpoint = 0.0;
#ifdef __linux__
			// Output in seconds
			return  ((REAL_EXT) ( trialEnd.tv_sec - trialStart.tv_sec )) + ((REAL_EXT) ( trialEnd.tv_nsec - trialStart.tv_nsec )) / ((REAL_EXT) 1e9);
#elif defined(__MACH__)
			// Compute time span
			trialEnd -= trialStart;

			// Output in seconds
			return ((REAL_EXT) trialEnd ) * tickRatio;
#elif defined(_WIN32)
			// Output in seconds
			return (trialEnd.QuadPart - trialStart.QuadPart) / tickRatio;
#endif
		}

		//! Output the status of the simulation
		void SimulationInfo::reportProgress(UNSIGNED_INTEGER done, UNSIGNED_INTEGER total)
		{
			if(NULL != ptrPSSA->ptrProgrCallback)
				(*(ptrPSSA->ptrProgrCallback))(done, total, ptrPSSA->ptrProgrCallbackUserData);
			else
			{
				ofsStatus.seekp(0);
				ofsStatus << "Progress ... " << done << '/' << total << "     " << std::endl << std::flush;
			}
		}

		//! Output the status of the simulation
		void SimulationInfo::reportError(STRING error)
		{
			// store error string
			strLastError = error;
			// try to report
			if((NULL != ptrPSSA)&&(NULL != ptrPSSA->ptrErrCallback))
				(*(ptrPSSA->ptrErrCallback))(error, ptrPSSA->ptrErrCallbackUserData);
			else
			{
				if(!bDisableOutput)
					std::cerr << error << std::endl;
				if(ofsStatus.is_open())
					ofsStatus << error << std::flush;
			}
			if(ofsDebug.is_open())
				ofsDebug << error << std::flush;
		}

		//! Provides output to file
		void SimulationInfo::doOutput()
		{
			// Check if output is disabled or scheduled to start later
			if((bDisableOutput)||(dTimeSimulation < dTimeCheckpoint))
				return;

			// Temporary variables
			REAL						dTemp;
			static std::stringstream	ssTemp;
			static STRING				sOutLine;

			// Calculate the time difference
			if(dTimeSimulation > dTimeEnd)
				dTemp = fabs(dTimeEnd - dTimeCheckpoint);
			else
				dTemp = fabs(dTimeSimulation - dTimeCheckpoint);

			if(dTemp >= dTimeStep)
			{
				// Output the population to stream
				for(UNSIGNED_INTEGER si = 0; si < ptrPSSA->ptrData->unSubvolumes; si++) {
					if (si != 0) ssTemp << "\t";
					for(UNSIGNED_INTEGER i = 1; i < ptrPSSA->ptrData->unSpecies; i++) {
						if (i != 1) ssTemp << ",";
						ssTemp << ptrPSSA->ptrData->arSubvolumes[si].aruN[i];
					}
				}
				ssTemp << std::endl;

				// Store to string & clear temporary stream
				sOutLine = ssTemp.str();
				ssTemp.str("");

				// Output to file
				do
				{
					// Advance in time
					dTimeCheckpoint	+= dTimeStep;
					dTemp			-= dTimeStep;

					// Output to file
					ssOutput << sOutLine;

					frame++;

					if (vtkOutName.size() != 0)
					{
						std::stringstream rank_str;
						rank_str << "_" << mpi.getRank() << "_";
						vtkProcess.WriteVtkFile(vtkOutName + rank_str.str() ,this,frame);
					}
					else
					{
						ssOutput << sOutLine << "\n";
					}
					
					reportProgress(dTimeCheckpoint * 100.0/dTimeEnd, 100.0);
				}
				while(dTemp >= dTimeStep);
			}
		}

		//! Called within the time-sampling routine to
		//! perform an update for a delayed reaction.
		//! FIXME : test non-consuming reactions
		bool SimulationInfo::UpdateCallback( void )
		{
			// Output the population
			doOutput();

			// Check reaction type
#if BT_TODO
			if(pssalib::datamodel::DataModel::D2 == ptrPSSA->ptrData->arReactClass[ptrPSSA->ptrData->mu])
				uUpdateMode = 2;	// products
			else if(pssalib::datamodel::DataModel::D1 == ptrPSSA->ptrData->arReactClass[ptrPSSA->ptrData->mu])
				uUpdateMode = 0;	// ALL
			else
			{
				reportError("SimulationInfo::UpdateCallback() - Error: unknown reaction type.");
				return false;
			}
#endif

			bOverrideUpdate = true;
			ptrPSSA->ptrUpdate->do_update(this);
			bOverrideUpdate = false;

			return true;
		}
	}
}
