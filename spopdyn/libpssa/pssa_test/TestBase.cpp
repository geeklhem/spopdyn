#include "TestBase.h"

void reaction_callback_wrapper(pssalib::datamodel::DataModel* dm, REAL t, void* user)
{
	TestBase* base = reinterpret_cast<TestBase*>(user);
	base->ReactionCallback(dm, t);
}

void progress_callback(UNSIGNED_INTEGER a, UNSIGNED_INTEGER b, void* user)
{
	std::cout << "\rProgress: " << a << "/" << b << "   " << std::flush;
}

TestBase::TestBase()
	: pssa(NULL)
{
}

TestBase::~TestBase()
{
	if (pssa) {
		delete pssa;
	}
}

bool TestBase::Test()
{
	if (!Setup()) {
		return false;
	}

	bool result = pssa->run_avg(&simInfo);
	std::cout << std::endl;
	return result;
}

bool TestBase::Setup()
{
	simInfo.dTimeStart = 0.0;
	simInfo.dTimeStep = 0.1;
	simInfo.dTimeEnd = 1000.0;
	simInfo.strOutput = "out/";
	simInfo.arSamples.push_back(10);
	simInfo.arSpeciesIds.push_back("all species");


	pssa = new pssalib::PSSA();
	if (!pssa) {
		std::cerr << "Failed to allocate the simulation engine." << std::endl;
		return false;
	}

	if(!pssa->setMethod(pssalib::PSSA::M_SPDM))	{
		std::cerr << "Failed to set simulation method." << std::endl;
		return false;
	}

	pssa->SetReactionCallback(&reaction_callback_wrapper, this);
	pssa->SetProgressCallback(&progress_callback, this);

	return true;
}

void TestBase::ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t)
{
}