#pragma once

#include "TestBase.h"

class TestDiffusion : public TestBase
{
public:
	TestDiffusion();
	virtual ~TestDiffusion();

private:
	virtual bool Setup();
	virtual void ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t);
};