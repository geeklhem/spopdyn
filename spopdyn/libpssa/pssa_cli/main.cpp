#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <sbml/SBMLTypes.h>
LIBSBML_CPP_NAMESPACE_USE

#include "PSSA.h"

namespace po = boost::program_options;

// Debug backtrace
#ifdef __linux__
#include <execinfo.h>
#include <signal.h>

void handler(int sig)
{
	void *array[10];
	size_t size;

	// get void*'s for all entries on the stack
	size = backtrace(array, 10);

	// print out all the frames to stderr
	fprintf(stderr, "Error: signal %d:\n", sig);
	backtrace_symbols_fd(array, size, 2);
	exit(1);
}
#endif

// Parse an array of numbers represented as strings
template< class T >
void parseArray(const STRING &sInput, std::vector<T> *arOut, const char* delim = ",", UNSIGNED_INTEGER nMax = 0) throw(boost::bad_lexical_cast &)
{
    std::vector<STRING>				vStrs;
    std::vector<STRING>::iterator	itS;

    // Limit number of elements
    UNSIGNED_INTEGER				n = 0;

    // Split the string & remove extra delimiters
    boost::split(vStrs,sInput,boost::is_any_of(delim));
    itS = std::remove(vStrs.begin(), vStrs.end(), STRING(""));
    vStrs.erase(itS, vStrs.end());

    for(itS = vStrs.begin(); (itS != vStrs.end()) && ((0 == nMax)||(n < nMax)); ++itS, ++n)
        arOut->push_back(boost::lexical_cast<T>((*itS)));
}

bool process_program_options(po::variables_map& vm, int argc, char** argv)
{
	try {
		po::options_description generic("Generic options");
		generic.add_options()
			("help,h",																"produce help message")
			("output_path,o",	po::value<STRING>(),								"Output path")
			("sbml_file,i",		po::value<STRING>(),								"SBML input file")
			("species,s",		po::value<STRING>(),								"Comma-separated list of species ids for which histograms will be computed")
			("dt",				po::value<REAL>()->default_value(0.1),				"Time interval between outputs")
			("tend",			po::value<REAL>()->default_value(1000.0),			"End time of the simulation")
			("tstart",			po::value<REAL>()->default_value(0.0),				"Time at which to begin outputing trajectories")
			("ntrajectories,n",	po::value<STRING>()->default_value("10"),			"Comma-separated list of number of trajectories to simulate")
			("mode",			po::value<STRING>()->default_value("stat"),			"Mode of the simulation, can be either \"trial\" or \"stat\"")
			("methods,m",		po::value<STRING>(),								"Comma-separated list of simulation method ids:" \
																					"\n0 - Gillespie's Direct Method" \
																					"\n1 - Partial Propensity Direct Method" \
																					"\n2 - PSSA with Composition-Rejection Sampling" \
																					"\n3 - Sorting Partial Propensity Direct Method")
			("stat_info",		po::value<STRING>(),					"In \"Stat\" mode, causes the engine to switch between one- and multidimensional histograms")
			("boundary,b",		po::value<STRING>()->default_value("periodic"),		"Boundary conditions, can be either \"periodic\" or \"reflexive\"")
			("volumetype",		po::value<STRING>()->default_value("single"),		"Volume description, can be either:" \
																					"\n\"single\" - a single subvolume (no diffusion)" \
																					"\n\"homogeneous2d\" - 2-dimensional homogeneous grid" \
																					"\n\"homogeneous3d\" - 3-dimensional homogeneous grid")
			("gridpoints",		po::value<UNSIGNED_INTEGER>()->default_value(4),	"Number of grid points along one dimension")
			("omega",			po::value<REAL>()->default_value(1.0),				"Size of the total volume")
			("initialpop",		po::value<STRING>()->default_value("distribute"),	"Specifies how the initial population specified in the SBML file is distributed across the volume:" \
																					"\n\"distribute\" - the population is evenly distributed, i.e. each subvolume gets pop/num_volumes" \
																					"\n\"concentrate\" - the population is concentrated in the middle cell, i.e. one subvolume gets everything" \
																					"\n\"multiply\" - the population is multiplied, i.e. each subvolume gets the total population")
			("verbose,v",									"\nIf specified, the program will print additional information to cout")
			("gnuplot",		po::value<STRING>()->default_value("yes"),							"\nvisualize result with gnuplot \"Stat\" mode")
			("vtk",			po::value<STRING>()->default_value("yes"),							"\nvisualize result with vtk \"Trial\" mode")
		;
	
		po::options_description cmdline("Command-line");
		cmdline.add_options()
			("config_file,c",	po::value<STRING>()->default_value("config.cfg"),	"Configuration file")
		;
		cmdline.add(generic);
	
		po::options_description cfgfile("Configuration");
		cfgfile.add(generic);

		// First parse the command line, to retrieve the config_file parameter.
		po::store(po::parse_command_line(argc, argv, cmdline), vm);

		// check for help options
	
		if (vm.count("help")) {
			std::cout << cmdline << std::endl;
			return false;
		}
		
		// Then parse the config file.
		std::string cfgfilename = vm["config_file"].as<std::string>();
		
		po::store(po::parse_config_file<char>(cfgfilename.c_str(), cfgfile), vm);

		// Now read the command line again, to override the config file parameters where applicable.
		po::store(po::parse_command_line(argc, argv, cmdline), vm);

		po::notify(vm);
		
	} catch (po::error &e) {
		std::cerr << "Error processing program options: " << e.what() << std::endl;
		return false;
	}

	return true;
}

int main(int argc, char** argv)
{
	pssalib::datamodel::SimulationInfo simInfo;
#ifdef __linux__
	signal(SIGSEGV, handler);
#endif
	if (simInfo.mpi.onlyMaster())
		std::cout << std::endl << "libpssa command line interface\n" << std::endl << std::endl;
	
	// Parse command line and the configuration file
	po::variables_map vm;
	if (!process_program_options(vm, argc, argv)) {
		return -1;
	}
	
	
	// Set up the simulation info and the pssa instance.
	pssalib::PSSA *pssa = new pssalib::PSSA();
	if (!pssa) {
		std::cerr << "Error : failed to allocate the simulation engine." << std::endl;
		return -2;
	}

	std::vector<UNSIGNED_INTEGER> methods;
	STRING inputFile;
	STRING outputPath;
	STRING mode;
	STRING stat_info;
	
	try {
		simInfo.bVerbose = vm.count("verbose") ? true : false;
		simInfo.dTimeStart = vm["tstart"].as<REAL>();
		simInfo.dTimeStep = vm["dt"].as<REAL>();
		simInfo.dTimeEnd = vm["tend"].as<REAL>();
		
		parseArray(vm["species"].as<STRING>(), &simInfo.arSpeciesIds);
		parseArray(vm["methods"].as<STRING>(), &methods);
		parseArray(vm["ntrajectories"].as<STRING>(), &simInfo.arSamples);
		inputFile = vm["sbml_file"].as<STRING>();
		outputPath = vm["output_path"].as<STRING>();
		mode = vm["mode"].as<STRING>();
		
		if (mode == "trial")
			simInfo.stStatistics = pssalib::datamodel::SimulationInfo::ST_SINGLE;
		else
			simInfo.stStatistics = pssalib::datamodel::SimulationInfo::ST_SINGLE;
		
		stat_info = vm["stat_info"].as<STRING>();
		
		simInfo.unNumGridPoints = vm["gridpoints"].as<UNSIGNED_INTEGER>();
		simInfo.dOmega = vm["omega"].as<REAL>();
		
		STRING strBoundary = vm["boundary"].as<STRING>();
		if (strBoundary == "periodic") simInfo.eBoundaryConditions = pssalib::datamodel::SimulationInfo::BC_Periodic;
		else if (strBoundary == "reflexive") simInfo.eBoundaryConditions = pssalib::datamodel::SimulationInfo::BC_Reflexive;

		STRING strVolume = vm["volumetype"].as<STRING>();
		if (strVolume == "single") simInfo.eVolumeType = pssalib::datamodel::SimulationInfo::VT_SingleSubvolume;
		else if (strVolume == "homogeneous2d") simInfo.eVolumeType = pssalib::datamodel::SimulationInfo::VT_Homogeneous2d;
		else if (strVolume == "homogeneous3d") simInfo.eVolumeType = pssalib::datamodel::SimulationInfo::VT_Homogeneous3d;

		STRING strInitialPop = vm["initialpop"].as<STRING>();
		if (strInitialPop == "distribute") simInfo.eInitialPop = pssalib::datamodel::SimulationInfo::IP_Distribute;
		else if (strInitialPop == "concentrate") simInfo.eInitialPop = pssalib::datamodel::SimulationInfo::IP_Concentrate;
		else if (strInitialPop == "multiply") simInfo.eInitialPop = pssalib::datamodel::SimulationInfo::IP_Multiply;
		
		if (vm["gnuplot"].as<STRING>() == STRING("yes"))
			simInfo.plotVisualize = true;
		else
			simInfo.plotVisualize = false;
		
		if (vm["vtk"].as<STRING>() == STRING("yes"))
			simInfo.vtkVisualize = true;
		else
			simInfo.vtkVisualize = false;
		
	} catch (po::error &e) {
		std::cerr << "Error processing program options: " << e.what() << std::endl;
		return -1; 
	}

	if (simInfo.arSamples.empty())
	{
		std::cerr << "Number of trajectories not specified" << std::endl;
		return -3;
	}
	if (methods.empty())
	{
		std::cerr << "No simulation method specified" << std::endl;
		return -3;
	}
	if (simInfo.arSpeciesIds.empty())
	{
		simInfo.arSpeciesIds.push_back("all species");
	}

	if (simInfo.bVerbose)
	{
		std::cout << "Simulate model defined in :" << std::endl << "\tinput file : " << inputFile << std::endl;
		std::cout << "and output results to" << std::endl << "\tsave path : " << outputPath << std::endl;
		std::cout << "using each of these methods:" << std::endl;
		for(std::vector<UNSIGNED_INTEGER>::iterator mi = methods.begin(); mi != methods.end(); mi++)
			std::cout << pssalib::PSSA::GetMethodName((pssalib::PSSA::EMethod)*mi) << '\t';
		std::cout << std::endl;
		
		if (mode == "trial")
		{
			std::cout << "Run " << simInfo.arSamples[0] << " trials until time = '" << simInfo.dTimeEnd << "'";
			std::cout << " outputting population after each '" << simInfo.dTimeStep << "' seconds beginning at t = '" << simInfo.dTimeStart << "'" << std::endl;
		}
		else if (mode == "stat")
		{
			std::cout << "Compute PDFs over each of following numbers of trajectories: " << std::endl;
			for(std::vector<UNSIGNED_INTEGER>::iterator ti = simInfo.arSamples.begin(); ti != simInfo.arSamples.end(); ti++)
				std::cout << (*ti) << '\t';
			std::cout << std::endl << "at time = '" << simInfo.dTimeEnd << "'" << std::endl;
		}

		if(simInfo.arSpeciesIds[0] != "all species")
		{
			std::cout << "Output results only for species with following ids :" << std::endl;
			for(std::vector<STRING>::iterator si = simInfo.arSpeciesIds.begin(); si != simInfo.arSpeciesIds.end(); ++si)
				std::cout << "'" << (*si) << "'\t";
			std::cout << std::endl;
		}
		else
		{
			std::cout << "Output results for all species" << std::endl;
		}

		std::cout << std::endl;
	}
	
	if (!simInfo.readSBMLFile(inputFile))
	{
		std::cerr << "Error : failed to load model file '" << inputFile << "'." << std::endl;
		delete pssa;
		return -2;
	}

	// Run all the simulations
	for(std::vector<UNSIGNED_INTEGER>::iterator mi = methods.begin(); mi != methods.end(); mi++)
	{
		if(!pssa->setMethod((pssalib::PSSA::EMethod)(*mi)))
		{
			std::cerr << "Error : failed to set simulation method " << *mi << std::endl;
			continue;
		}
		
		// Create the output directory
		STRING path;
		{
			std::stringstream ssTemp;
			ssTemp << outputPath << '/' << pssalib::PSSA::GetMethodName((pssalib::PSSA::EMethod)*mi) << '/';
			path = ssTemp.str();
		}
		if (!pssalib::PSSA::MakeDir(path))
		{
			std::cerr << "Error : failed to create output directory" << std::endl;
			break;
		}
		simInfo.strOutput = path;

		if (simInfo.mpi.onlyMaster())
			std::cout << "simulating '" << pssa->getModelName() << "' using " << pssalib::PSSA::GetMethodName((pssalib::PSSA::EMethod)*mi) << "  ... " << std::endl << std::flush;
		
		// Decide which function to use
		bool result;
		if (mode == "trial")
		{
			simInfo.stStatistics = pssalib::datamodel::SimulationInfo::ST_SINGLE;
			result = pssa->run_avg(&simInfo);
		}
		else if (mode == "stat")
		{
			
			if(stat_info == STRING("multi"))	// multidimensional?
				simInfo.stStatistics = pssalib::datamodel::SimulationInfo::ST_MULTI;
			else
				simInfo.stStatistics = pssalib::datamodel::SimulationInfo::ST_SINGLE;
		  
			result = pssa->run_hist(&simInfo);
		}
		else
		{
			std::cerr << "Error : unknown mode " << mode << std::endl;
			break;
		}

		if (simInfo.mpi.onlyMaster())
		{
			if(result)
				std::cout << "done!" << std::endl;
			else
				std::cout << "FAILED!\n" << pssa->getLastError() << std::endl;
		}
	}

	delete pssa;

#if defined(_DEBUG) && defined(WIN32)
	std::cin.get();
#endif

	return 0;
}
