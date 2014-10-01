/*
 * 
 *    Gray Scott 2D, 
 * 
 *       Partial Propensity example
 * 
 *    Author: Pietro Incardona
 * 
 *    Warning this is an example code, results an performance are not related with any MOSAIC paper
 * 
 * */
#include "PSSA.h"
#include <sbml/SBMLTypes.h>
LIBSBML_CPP_NAMESPACE_USE

void init_callback(pssalib::datamodel::DataModel::SubvolumesInit * sb);
void reaction_callback_wrapper(pssalib::datamodel::DataModel* dm, REAL t, void* user);
void progress_callback(UNSIGNED_INTEGER a, UNSIGNED_INTEGER b, void* user);

#define U	1e6
#define H	0.01
#define N	64

#define UHH	U*H*H

pssalib::datamodel::SimulationInfo simInfo;

int main(int argc, char** argv)
{
	std::cout << "Running Gray Scott 2D example ... " << std::endl;
	
	std::string inputFile = "sbml/grayscott_sp.sbml";
	if (!simInfo.readSBMLFile(inputFile))
	{
		std::cerr << "Failed to load model file '" << inputFile << "'." << std::endl;
		return false;
	}

	// setup seed
	
	srand((unsigned)time(NULL));
	
	// Boundary condition periodic
	
	simInfo.eBoundaryConditions = pssalib::datamodel::SimulationInfo::BC_Periodic;
	
	// 2D Homogeneous space
	
	simInfo.eVolumeType = pssalib::datamodel::SimulationInfo::VT_Homogeneous2d;
	
	// Volume
	
	simInfo.dOmega = H*H*N*N;
	
	// number of points in the grid (one side)
	
	simInfo.unNumGridPoints = N;
	
	// User defined concentration, and call-back
	
	simInfo.eInitialPop = pssalib::datamodel::SimulationInfo::IP_UserDefined;
	simInfo.eInitialUserDefined = init_callback;
	
	// Starting time
	
	simInfo.dTimeStart = 0.0;
	
	// Time step between two configuration save
	
	simInfo.dTimeStep = 2000/U/U/1500;
	
	// End time of simulation
	
	simInfo.dTimeEnd = 2000/U/U*1.5;
	
	// output location
	
	simInfo.strOutput = "out";
	simInfo.arSamples.push_back(1);
	
	// out all species
	
	simInfo.arSpeciesIds.push_back("all species");
	
	// output name sequence
	
	simInfo.vtkOutName = "out/gray";
	
	// run VTK at the end od the simulation
	
	simInfo.vtkVisualize = false;

	
	// Create an istance of the simulation class
	
	pssalib::PSSA *pssa = new pssalib::PSSA();
	if (!pssa)
	{
		std::cerr << "Failed to allocate the simulation engine." << std::endl;
		return false;
	}

	// Set up the method to use to SPDM: Sort Partial Propensity Direct Method
	
	if(!pssa->setMethod(pssalib::PSSA::M_SPDM))
	{
		std::cerr << "Failed to set simulation method." << std::endl;
		return false;
	}

	// Set a callback for each reaction event
	
	pssa->SetReactionCallback(&reaction_callback_wrapper, NULL);
	
	// Set a callback for each 
	
	pssa->SetProgressCallback(&progress_callback, NULL);

	// Run the simulation
	
	bool result = pssa->run_avg(&simInfo);
	
	// end
	
	std::cout << std::endl;
	
	return 0;
}

void init_callback(pssalib::datamodel::DataModel::SubvolumesInit * sb)
{
	int start_sp1 = 0;
	int start_sp2 = 0;
	
	// Set initial population
  
	// Gray Scott simulation require two species + 1 (source/sink)
	
	if ( simInfo.getNumSpecies() != 2 + 1 )
		throw "init_callback: Error Gray Scott example work with only two species";
	
	// For each subvolume
	
	for (int i = 0 ; i < sb->getNumSub() ; i++)
	{
		// random var from 0 to 1
	  
		double r = ((double)rand())/RAND_MAX;
	  
		int s1 = i % 64;
		int s2 = i/64;
		
		if (s1 > 24 && s1 < 40 && s2 > 24 && s2 < 40)
		{start_sp1 = UHH/2 + (0.04*(r-0.5)*UHH + 0.5);
		 start_sp2 = UHH/4 + (0.02*(r-0.5)*UHH + 0.5);}
		else
		{start_sp1 = UHH;
		 start_sp2 = 0;}
		
		// Set population for species 1,2 subvolumes i

		sb->setInit(i,1,start_sp1);
		sb->setInit(i,2,start_sp2);
	}
}

void reaction_callback_wrapper(pssalib::datamodel::DataModel* dm, REAL t, void* user)
{
	// do nothing
}

int total_frame = 0;

void progress_callback(UNSIGNED_INTEGER a, UNSIGNED_INTEGER b, void* user)
{
	// display progress
  
	std::cout << "\rProgress: " << a << "/" << b << " frame ( " << total_frame << " )" << std::flush;
	total_frame++;
}

