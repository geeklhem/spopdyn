#pragma once

#include "TestBase.h"

class TestReaction : public TestBase
{
public:
	TestReaction();
	virtual ~TestReaction();

private:
	virtual bool Setup();
	virtual void ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t);
};