#include "TestReactionDiffusion.h"

TestReactionDiffusion::TestReactionDiffusion()
{
}

TestReactionDiffusion::~TestReactionDiffusion()
{
}

bool TestReactionDiffusion::Setup()
{
	std::string inputFile = "sbml/Multimerization.sbml";
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

void TestReactionDiffusion::ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t)
{
	INTEGER total = 0;
	INTEGER max = 0;
	for (UNSIGNED_INTEGER i = 0; i < dm->unSubvolumes; i++) {
		pssalib::datamodel::DataModel::Subvolume& sv = dm->arSubvolumes[i];
		INTEGER thiscount = 0;
		thiscount += sv.aruN[1];
		thiscount += sv.aruN[2] * 2;
		thiscount += sv.aruN[3] * 3;
		thiscount += sv.aruN[4] * 4;
		thiscount += sv.aruN[5] * 5;
		total += thiscount;
		if (thiscount > max) max = thiscount;
	}

	// Totals shouldn't change.
	if (total != 100) {
		std::cerr << "TestReactionDiffusion::ReactionCallback: Incorrect number of monomers, total=" << total << std::endl;
	}

	// By t=1, molecules should have distributed among the subvolumes.
	if (t > 1.0) {
		if (max > 70 || max > 70) {
			std::cerr << "TestReactionDiffusion::ReactionCallback: Molecules do not seem to be diffusing, max=" << max << std::endl;
		}
	}

	TestBase::ReactionCallback(dm, t);
}