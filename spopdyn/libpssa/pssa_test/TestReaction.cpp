#include "TestReaction.h"

TestReaction::TestReaction()
{
}

TestReaction::~TestReaction()
{
}

bool TestReaction::Setup()
{
	std::string inputFile = "sbml/Multimerization.sbml";
	if (!simInfo.readSBMLFile(inputFile))
	{
		std::cerr << "Failed to load model file '" << inputFile << "'." << std::endl;
		return false;
	}

	simInfo.eBoundaryConditions = pssalib::datamodel::SimulationInfo::BC_Periodic;
	simInfo.eVolumeType = pssalib::datamodel::SimulationInfo::VT_SingleSubvolume;
	simInfo.dOmega = 1.0;
	simInfo.unNumGridPoints = 1;
	simInfo.eInitialPop = pssalib::datamodel::SimulationInfo::IP_Concentrate;

	return TestBase::Setup();
}

void TestReaction::ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t)
{
	INTEGER total = 0;
	for (UNSIGNED_INTEGER i = 0; i < dm->unSubvolumes; i++) {
		pssalib::datamodel::DataModel::Subvolume& sv = dm->arSubvolumes[i];
		total += sv.aruN[1];
		total += sv.aruN[2] * 2;
		total += sv.aruN[3] * 3;
		total += sv.aruN[4] * 4;
		total += sv.aruN[5] * 5;
	}

	// Totals shouldn't change.
	if (total != 100) {
		std::cerr << "TestReaction::ReactionCallback: Incorrect number of monomers, total=" << total << std::endl;
	}

	TestBase::ReactionCallback(dm, t);
}