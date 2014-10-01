#pragma once

#include "PSSA.h"

class TestBase
{
public:
	TestBase();
	virtual ~TestBase();

	bool Test();

protected:
	virtual bool Setup();
	virtual void ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t);

	pssalib::datamodel::SimulationInfo simInfo;
	pssalib::PSSA *pssa;

	friend void reaction_callback_wrapper(pssalib::datamodel::DataModel* dm, REAL t, void* user);
};