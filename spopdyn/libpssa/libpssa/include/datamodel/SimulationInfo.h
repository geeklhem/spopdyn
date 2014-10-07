/**
 * \file SimulationInfo.h
 * \brief Declares a datatype that is used to marshall simulation data between library modules.
 *
 * \date   03.06.2011
 * \author sanyi
 */

#ifndef SIMULATIONINFO_H_
#define SIMULATIONINFO_H_

// Forward declarations
namespace pssalib
{
	class PSSA;
	namespace datamodel
	{
		class DataModel;
	}
}

#include "./../PSSA.h"
#include "./../vtk/vtk.h"
#include "./../gnuplot/IBgnuplot.h"
#include "./../MPI/mpi_spread.h"

namespace pssalib
{
	namespace datamodel
	{	  
		////////////////////////////////////////
		//! \brief Simulation information class
		class SimulationInfo
		{
		/////////////////////////////////
		// Constructors
		public:
			// Constructor
			SimulationInfo();

			// Copy constructor
			SimulationInfo(SimulationInfo & right);

			// Destructor
			~SimulationInfo();

		/////////////////////////////////
		// Typedefs
		public:
			// Data type for statistic output
			typedef enum tagStatisticsType
			{
				ST_NONE,	//!<No statistical output
				ST_SINGLE,	//!<"Trial" mode: one-dimensional histogram.
							//! "Stat" mode: Statistics (std. deviation & mean) for each species separately
				ST_MULTI	//!<"Stat" mode: Output data for a multidimensional histogram
			} StatisticsType;

			typedef enum tagVolumeType
			{
				VT_Invalid,		//!<Invalid Volume Type
				VT_SingleSubvolume,	//!<Only one subvolume
				VT_Homogeneous2d,	//!<Homogeneous 4-Way connectivity
				VT_Homogeneous3d	//!<Homogeneous 6-Way connectivity 
			} VolumeType;

			typedef enum tagBoundaryConditionsType
			{
				BC_Invalid,		//!<Invalid Boundary condition
				BC_Reflexive,		//!<Reflexive Boundary condition
				BC_Periodic		//! Periodic Boundary condition
			} BoundaryConditionsType;

			typedef enum tagInitialPopulationType
			{
				IP_Invalid,		//!<Invalid initial population
				IP_Distribute,		//!<Uniformly distributed initial population
				IP_Concentrate,		//!<Population concentrated in the middle
				IP_Multiply,		//!<Get the initial amount from the SBML file and uses this value for all cells
				IP_UserDefined,		//!<User Defined (need callback)
				IP_Sbml //!<get the values from the SBML annotation
			} InitialPopulationType;
			
		/////////////////////////////////
		// Attributes
		protected:
			//! Pointer to simulation engine class.
			pssalib::PSSA*			ptrPSSA;

#ifdef __linux__
			struct timespec trialStart, trialEnd;	//!< Internal timing variables.
#elif defined(__MACH__)
			UINT64			trialStart, trialEnd;	//!< Internal timing variables.
			REAL_EXT		tickRatio;				//!< Internal timing variables.
#elif defined(_WIN32)
			LARGE_INTEGER	trialStart, trialEnd;	//!< Internal timing variables.
			REAL_EXT		tickRatio;				//!< Internal timing variables.
#endif

			//! True if update mode is overridden by the callback.
			bool					bOverrideUpdate;

			//! Last error string
			STRING					strLastError;
			
			//! frame number
			int 					frame;
			
		public:

			//! MPI [RESERVED]
		  
			pssalib::mpi::mpi_spread mpi;
		  
			//! Gnuplot [RESERVED]
		  
			pssalib::gnuplot::IBgnuplot plot;
			
			//! Visualize or not with gnuplot the output after the simulation [IN OPTIONAL, dafault: false]
			bool plotVisualize;
		  
			//! vtkProcess [RESERVED]
		  
			pssalib::vtk::vtk vtkProcess;
		  
			//! Pointer to a model [RESERVED]
			SBMLDocument	*ptrSBMLDocument;

			//! Output path (directory) [IN MADATORY]
			STRING					strOutput;
			
			//! vtk out filename [IN MADATORY]
			STRING					vtkOutName;

			//! Visualize or not with paraview the output after the simulation [IN OPTIONAL, dafault: false]
			bool vtkVisualize;
			
			//! Array of species to output for the probability density [IN OPTIONAL, default: "all species"]
			std::vector<STRING>		arSpeciesIds;

			//! Array of number of samples for each simulation [IN OPTIONAL, default arSamples[0] = 1]
			std::vector<UNSIGNED_INTEGER>		arSamples;

			// Simulation timing
			REAL						dTimeCheckpoint,		//!<last output time [RESERVED]
									dTimeStart,			//!<initial output time [IN OPTIONAL, default = 0.0]
									dTimeStep,			//!<time interval between subsequent outputs [IN MANDATORY]
									dTimeEnd,			//!<end time of the simulation [IN MANDATORY]
									dTimeSimulation;		//!<current simulation time [RESERVED]

			//! Type of statistic output [IN OPTIONAL, default ST_NONE]
			StatisticsType			stStatistics;

			//! Geometry description, see tagVolumeType [IN MADATORY]
			VolumeType				eVolumeType;
			
			//! Boundary description, see tagBoundaryConditionsType [IN MADATORY]
			BoundaryConditionsType	eBoundaryConditions;
			
			//! Size of one side of the grid (if unNumGridPoints = N; 2D = N*N 3D = N*N*N) [IN MADATORY]
			UNSIGNED_INTEGER		unNumGridPoints;
			
			//! Volume of one cell [IN MADATORY]
			REAL					dOmega;
			
			//! Initial Population Type see tagInitialPopulationType [IN MADATORY]
			InitialPopulationType	eInitialPop;
			
			//! Initial population callback [IN Otional, dafault: NULL]
			void (* eInitialUserDefined)(DataModel::SubvolumesInit * sub);

			//! Update mode [RESERVED]
			// 0 - Update ALL
			// 1 - Update reactants
			// 2 - Update products
			UNSIGNED_INTEGER		uUpdateMode;

			//! Output [RESERVED]
			std::stringstream		ssOutput;
			//! Error [RESERVED]
			std::stringstream		ssError;
			//! Debug [RESERVED]
			std::ofstream			ofsDebug;
			//! Status [RESERVED]
			std::ofstream			ofsStatus;

			// Status file name

			// flags : disable output, verbose output, statistics
			// statistics - in "Trial" mode, computes the additional statistics (see help)
			//              in "Stat" mode, controls the output
			
			//! Disable output [IN OPTIONAL: default no]
			bool					bDisableOutput;
			
			//! Verbosity [IN OPTIONAL: default no]
			bool					bVerbose;

			//////////////////////////////
			// Methods
			protected:
				//! \internal Check if given path is valid
				bool	chkPath(STRING& path) const
				{
					std::fstream	fTest;

		#if defined(__linux__) || defined(__MACH__)
					if(path[path.length()-1] != '/')
						path += '/';
		#elif defined(_WIN32)
					std::replace(path.begin(), path.end(), '\\', '/');
					if(path[path.length()-1] != '/')
						path += '/';
		#else
			#error	Write a suitable version of code above for your platform.
		#endif
					STRING	strDate, strFilePath;
					getDateStr(strDate);
					strFilePath = path + strDate;

					fTest.open(strFilePath.c_str(), std::fstream::out);

					if(fTest.is_open())
					{
						fTest.close();
						remove(strFilePath.c_str());
						return true;
					}
					else
					{
						return false;
					}
				};

				//! \internal Returns current date and time (used for file/folder names)
				void	getDateStr(STRING& str) const
				{
					static char			buf[128];
					static struct tm	*today;

					time_t calender_time = time(NULL);
					today = localtime( &calender_time );

					strftime ( buf, 128, "%d-%m-%y_%H-%M-%S", today );

					str = buf;
				};

		public:

			/////////////////////////////////
			// Simulation engine interface

			// Attach to simulation engine
inline		void			attachPSSA(PSSA*	pPSSA)	{	ptrPSSA = pPSSA;	};

			// Detach to simulation engine
inline		void			detachPSSA()				{	ptrPSSA = NULL;		};

			// Check if this instance is attache to a simulation engine
inline		bool			isAttached()		const	{	return (bool)(NULL != ptrPSSA);	};

			// Check if this instance is valid
			bool			isValid();

			// Check if data was loaded (or Test function is set)
inline		bool			isDataLoaded()		const	{	return	(bool)((NULL != ptrSBMLDocument)&&(ptrSBMLDocument->getNumErrors() == 0));	};


			/////////////////////////////////
			// Getters & setters

			//! Return number of species
inline 		int 	getNumSpecies()			{if (getDataModel() != 0)	return getDataModel()->unSpecies;
							 else {throw "SimulationInfo::getNumSpecies(): Error DataModel not loaded";}}

			//! Check whether the update mode is set externally
inline		bool			getOverrideUpdate()	const	{	return bOverrideUpdate;	};

			//! Get the underlying PSSA object
inline		PSSA*			getPSSA()			const	{	return ptrPSSA;	};

			//! Get the data model associated with PSSA object
inline		DataModel*		getDataModel()		const	{	return ptrPSSA->ptrData;	};

			//! Returns last error string
inline		STRING			getLatError()		const	{	return strLastError;	};

			//! Get the SBML model
inline		Model*	getSBMLModel()		const	{	return ptrSBMLDocument->getModel();	};

			/////////////////////////////////
			// Functions used within the
			// simulation loop

			//! True if the simulation has not reached dTimeEnd
inline		bool			isRunning()			const	{	return (dTimeSimulation < dTimeEnd);	};

			// Output the results to file
			void			doOutput( void );

			// FIXME : callback to update delayed reactions (consuming AND non-consuming(!) )
			bool			UpdateCallback( void );

			// load a new model file
			bool			readSBMLFile( STRING strFile );

			/////////////////////////////////
			// Report status

			// Output the status of the simulation
			void			reportProgress(UNSIGNED_INTEGER done, UNSIGNED_INTEGER total);

			// Output the status of the simulation
			void			reportError(STRING strError);

			/////////////////////////////////
			// Functions called before and
			// after each trial for timing

			// Initialize the timing
			void			beginTrial();

			// Calculate time spent on the run
			// Returns the time spent of this trial
			REAL_EXT		endTrial();
		};
	}
}

#endif
