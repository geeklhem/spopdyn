/**
 * \file mpi_spread.h
 * \brief small class to spread processes along multiple indipendent trajectories and collect data result. Input data consistency between process is outside the porpouse of this class
 *
 * \date   24.2.2013
 * \author Pietro Incardona
 */

#include "config.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

 namespace pssalib
{
	namespace mpi
	{
	  
		class mpi_spread
		{
			private:

				bool start_it;
				int rank;
				int n_proc;
				long int block;

			public:
			  
				mpi_spread();
				~mpi_spread();
			  
				void closeOut();
				void redirectOut();
				
				int getRank();
				
				unsigned int dataChuckSize(unsigned int size);
				
				void init(long int size);
				
				bool spread(unsigned int & tr);
				
				void collect_result(void * sbuf, void * rbuf, long int num_res, long int size_res);
				
				bool onlyMaster();
		};
	  
	}
}