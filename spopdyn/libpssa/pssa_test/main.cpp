#include <sbml/SBMLTypes.h>
LIBSBML_CPP_NAMESPACE_USE


#include "TestDiffusion.h"
#include "TestReaction.h"
#include "TestReactionDiffusion.h"

int main(int argc, char** argv)
{
	std::cout << "Running pure diffusion test..." << std::endl;

	TestBase* test_diffusion = new TestDiffusion;
	if (!test_diffusion->Test()) {
		std::cerr << "Diffusion test failed!" << std::endl;
	}
	
	std::cout << "Running pure reaction test..." << std::endl;
	
	TestBase* test_reaction = new TestReaction;
	if (!test_reaction->Test()) {
		std::cerr << "Reaction test failed!" << std::endl;
	}
	
	
	std::cout << "Running reaction-diffusion test..." << std::endl;

	
	TestBase* test_reactiondiffusion = new TestReactionDiffusion;
	if (!test_reactiondiffusion->Test()) {
		std::cerr << "Reaction-diffusion test failed!" << std::endl;
	}


	return 0;
}
