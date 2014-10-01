#pragma once

#include "TestBase.h"

class TestReactionDiffusion : public TestBase
{
public:
	TestReactionDiffusion();
	virtual ~TestReactionDiffusion();

private:
	virtual bool Setup();
	virtual void ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t);
};