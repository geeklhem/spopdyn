//IBS

#ifndef IBGNUPLOT_H_
#define IBGNUPLOT_H_

#include "config.h"


#include <fstream>
#include <iostream>
#include <string.h>
#include <cstdlib>
#include <vector>


namespace pssalib
{
	namespace gnuplot
	{
	  
		struct hist
		{
			int x;
			REAL p;
		};
	  
		class IBgnuplot
		{
			private:

			public:

				void GraphHist(STRING name);
		
				IBgnuplot();
				~IBgnuplot();
		};
	}
};

#endif
