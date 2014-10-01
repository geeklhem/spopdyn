#include "../../include/MPI/mpi_spread.h"
#include <fstream>

 namespace pssalib
{
	namespace mpi
	{
 
		std::ofstream dev_null;
		std::streambuf * backup_out;
		std::streambuf * backup_err;
	  
		// getRank
		
		int mpi_spread::getRank()
		{
			return rank;
		}
		
		// MPI close to dev/null
	  
		void mpi_spread::redirectOut()
		{
			#ifdef HAVE_MPI
		  
			if (rank != 0)
			{
				std::streambuf *psbuf;
				dev_null.open("/dev/null");

				backup_out = std::cout.rdbuf();     // back up cout's streambuf
				backup_err = std::cout.rdbuf();     // back up cout's streambuf

				psbuf = dev_null.rdbuf();   // get file's streambuf
				std::cout.rdbuf(psbuf);         // assign streambuf to cout
				std::cerr.rdbuf(psbuf);         // assign streambuf to cout
				
			}
			
			#endif
		}
		
		// MPI restore cout
		
		void mpi_spread::closeOut()
		{
			#ifdef HAVE_MPI
		  
			std::cout.rdbuf(backup_out);        // restore cout's original streambuf
			std::cerr.rdbuf(backup_err);        // restore cerr's original streambuf
			
			dev_null.close();
			
			#endif
		}
		
	  
		// Initialize MPI
	  
		mpi_spread::mpi_spread()
		:start_it(false)
		{
			n_proc = 1;
			rank = 0;
			
			#ifdef HAVE_MPI
		  
			MPI_Init(NULL,NULL);
			
			MPI_Comm_size(MPI_COMM_WORLD,&n_proc);
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		  
			#endif
		}
	  
		// Finalize
	  
		mpi_spread::~mpi_spread()
		{
			#ifdef HAVE_MPI
			
			MPI_Finalize();
			
			#endif
		}
	  
		//! Distribute work to each process
	  
		void mpi_spread::init(long int size)
		{
			start_it = false;
		  
			#ifdef HAVE_MPI
		  
			block = (int)size/n_proc;
			int thr = size%n_proc;
			
			// Increase the block, for the rest of trajectory
			
			if (rank < thr)
				block++;
			
			#else

			block = size;
				
			#endif
		}

		//! Return the operation chunk size assigned by the function spread() for each processor
		
		unsigned int mpi_spread::dataChuckSize(unsigned int size)
		{
			int block = (int)size/n_proc;
			int thr = size%n_proc;
			
			// Increase the block, for the rest of trajectory
			
			if (rank < thr)
				block++;
		  
			return block;
		}
		
		//! Spread like a parallel for, tr is the cicle counter, return true (still work to do) or false (no work to do)
		//! size of the cicle is initialized with init().
		
		bool mpi_spread::spread(unsigned int & tr)
		{
			if (start_it == false)
			{tr = rank*block;
			start_it = true;}
			else
				tr++;;
			  
			if (tr >= block + rank*block)
			{
				start_it = false;
				return false;
			}
			else
				return true;
		}

		void mpi_spread::collect_result(void * sbuf, void * rbuf, long int num_res, long int size_res)
		{
			// block can be different we cannot use Gatther
		  
			#ifdef HAVE_MPI

			if (sbuf == NULL || rbuf == NULL)	return;
			
			if (rank == 0)
			{
				unsigned int offs = block*size_res;
				for (int i = 1 ; i < n_proc ; i++)
				{
					int r_block = (int)num_res/n_proc;
					int thr = num_res%n_proc;
			
					// Increase the block, for the rest of trajectory

					if (i < thr)
						r_block++;
					
					MPI_Recv(&((char *)rbuf)[offs],r_block*size_res,MPI_CHAR,i,0,MPI_COMM_WORLD,NULL);
					offs += r_block*size_res;
				}
			}
			else
			{
				MPI_Send(sbuf,block*size_res,MPI_CHAR,0,0,MPI_COMM_WORLD);
			}
			
			#endif
		}

		bool mpi_spread::onlyMaster()
		{
			if (rank == 0)
				return true;
			else
				return false;
		}
		
	}
}