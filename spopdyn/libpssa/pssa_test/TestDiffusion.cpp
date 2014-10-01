#include "TestDiffusion.h"

TestDiffusion::TestDiffusion()
{
}

TestDiffusion::~TestDiffusion()
{
}

bool TestDiffusion::Setup()
{
	std::string inputFile = "sbml/Diffusion.sbml";
	if (!simInfo.readSBMLFile(inputFile))
	{
		std::cerr << "Failed to load model file '" << inputFile << "'." << std::endl;
		return false;
	}

	simInfo.eBoundaryConditions = pssalib::datamodel::SimulationInfo::BC_Periodic;
	simInfo.eVolumeType = pssalib::datamodel::SimulationInfo::VT_Homogeneous2d;
	simInfo.dOmega = 1.0;
	simInfo.unNumGridPoints = 3;
	simInfo.eInitialPop = pssalib::datamodel::SimulationInfo::IP_Concentrate;

	return TestBase::Setup();
}

void TestDiffusion::ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t)
{
	INTEGER totalA = 0, totalB = 0;
	INTEGER maxA = 0, maxB = 0;
	for (UNSIGNED_INTEGER i = 0; i < dm->unSubvolumes; i++) {
		pssalib::datamodel::DataModel::Subvolume& sv = dm->arSubvolumes[i];
		totalA += sv.aruN[1];
		totalB += sv.aruN[2];
		if (sv.aruN[1] > maxA) maxA = sv.aruN[1];
		if (sv.aruN[2] > maxB) maxB = sv.aruN[2];
	}

	// Totals shouldn't change.
	if (totalA != 100 || totalB != 100) {
		std::cerr << "TestDiffusion::ReactionCallback: Incorrect number of molecules, A=" << totalA << " B=" << totalB << std::endl;
	}

	// By t=1, molecules should have distributed among the subvolumes.
	if (t > 1.0) {
		if (maxA > 50 || maxB > 50) {
			std::cerr << "TestDiffusion::ReactionCallback: Molecules do not seem to be diffusing, maxA=" << maxA << " maxB=" << maxB << std::endl;
		}
	}

	TestBase::ReactionCallback(dm, t);
}