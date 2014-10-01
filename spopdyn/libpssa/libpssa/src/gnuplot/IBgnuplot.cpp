//IBS

#include "config.h"
#include "../../include/typedefs.h"
#include "../../include/gnuplot/IBgnuplot.h"

namespace pssalib
{
	namespace gnuplot
	{
		IBgnuplot::IBgnuplot()
		{
		}

		IBgnuplot::~IBgnuplot()
		{
		}
		
		void IBgnuplot::GraphHist(STRING name)
		{
			STRING plot_command;
			plot_command = STRING("gnuplot -e 'set boxwidth 0.5; set style fill solid 1.0 noborder; plot \"" + name + "\" using 1:2 title \"Probability\"  with boxes ; pause(-1);'");
			std::cout << "\n Plotting command:, " << plot_command << "\n";
			std::cout << "Plotting, press return to close \n";
			system(plot_command.c_str());
		}
	}
}

