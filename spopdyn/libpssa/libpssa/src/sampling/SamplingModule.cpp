/*
 * SamplingModule.cpp
 *
 *  Created on: 14 груд. 2010
 *      Author: sanyi
 */

#include "config.h"
#include "../../include/stdheaders.h"
#include "../../include/sampling/SamplingModule.h"
#include "../../include/datamodel/DataModel.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../dependencies/dSFMT-src-2.1/gsl_dSFMT.h"
#include <sstream>
#include "sha1.h"

namespace pssalib
{
	namespace sampling
	{
		SamplingModule::SamplingModule()
		{
			int rank = 1;
			std::stringstream pipe;
			
			#ifdef HAVE_MPI
			
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			rank++;
			
			#endif
			
			long int test = time(NULL)*getpid()*rank;
			
			SHA1* sha1 = new SHA1();
			sha1->addBytes( (const char *)&test, sizeof( long int ) );
			pipe << *(long int *)sha1->getDigest();
			pipe >> test;
			
			gsl_rng_default_seed = test/*(unsigned long)time(NULL)*/;
			ptrRNG = gsl_rng_alloc(gsl_rng_default);
			set_rng_seed(gsl_rng_default_seed);
		}

		SamplingModule::~SamplingModule()
		{
			gsl_rng_free(ptrRNG);
		}

		void SamplingModule::set_rng_seed(UNSIGNED_INTEGER seed)
		{
			gsl_rng_set(ptrRNG, seed);
		}

		bool SamplingModule::get_sample(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Sample time
			if(!sample_time(ptrSimInfo))
				return false;

			// Sample volume
			if(!sample_volume(ptrSimInfo))
				return false;

			// Sample reaction
			if(!sample_reaction(ptrSimInfo))
				return false;

			// Sample diffusion destination if necessary.
			pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();
			if(ptrData->IsDiffusionReaction(ptrData->mu))
			{
				// With isotropic diffusion we don't need a linear search.
				pssalib::datamodel::DataModel::Subvolume &sv = ptrData->arSubvolumes[ptrData->sv];
				UNSIGNED_INTEGER dstsv = (UNSIGNED_INTEGER)(gsl_rng_uniform_pos(ptrRNG) * sv.unDiffusionReactions);
				ptrData->svDst = sv.arNeighborVolumes[dstsv];
			}

			return true;
		}

		bool SamplingModule::sample_time(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();

			// Generate a random number
			
			REAL r;
			
			do
			{
				r = gsl_rng_uniform (ptrRNG);
			} while (r == 0);

			if(ptrData->vQueuedReactions.empty())
			{
				// Sample the time
				if(ptrData->dTotalPropensity <= 0.0)
				{
					// BT: uh; r == 0.0 does not mean we reached an absorbing state, right?
					// isn't it better to resample r; or use log1p(-r) instead of log(r)?
					ptrSimInfo->dTimeSimulation = std::numeric_limits<REAL>::infinity();
					return false;	// We have reached an absorbing state - exit
				}
				else
					ptrSimInfo->dTimeSimulation -= log(r) / ptrData->dTotalPropensity;

				return true;
			}
			else
			{
				REAL	T1, T2, at, F;
				std::vector<pssalib::datamodel::DataModel::DelayedReaction>::iterator curReaction = ptrData->vQueuedReactions.begin();

				// Initialize
				T1 =	ptrSimInfo->dTimeSimulation;
				T2 =	curReaction->time;
				at =	ptrData->dTotalPropensity * (T2 - T1);
				F  =	1.0 - exp(-at);

				while(F < r)
				{
					// Reached the end of the simulation
					ptrData->mu = curReaction->index;
					ptrSimInfo->dTimeSimulation = T2;
					// Write to file & update
					if(!ptrSimInfo->UpdateCallback())
						return false;

					// Check if we're still within simulation timespan
					if(T2 > ptrSimInfo->dTimeEnd)
						return false;

					// Remove the delayed reaction we just fired
					ptrData->vQueuedReactions.erase(curReaction);

					T1 = T2;
					// If all delayed reaction have already been fired
					// just calculate the time-step and exit.
					if(ptrData->vQueuedReactions.empty())
					{
						ptrSimInfo->dTimeSimulation = T1 - (boost::math::log1p(-r) + at) / ptrData->dTotalPropensity;
						return true;
					}

					curReaction = ptrData->vQueuedReactions.begin();
					T2 =	curReaction->time;
					at +=	ptrData->dTotalPropensity * (T2 - T1);
				    F  =	1.0 - exp(-at);
				}

				ptrSimInfo->dTimeSimulation = T2 - (boost::math::log1p(-r) + at) / ptrData->dTotalPropensity;
				return true;
			}

			return false;
		}

		bool SamplingModule::sample_volume(pssalib::datamodel::SimulationInfo* ptrSimInfo)
		{
			// Cast the data model to a suitable type
			pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();

			bool success;
			if (ptrData->unSubvolumes == 1) {
				ptrData->sv = 0;
				success = true;
			} else {
				UNSIGNED_INTEGER i;
				REAL r;
				success = crVolumeSampler.Sample(&ptrData->crsdVolume, ptrRNG, ptrData->dTotalPropensity, i, r);
				if (!success) {
					ptrData->sv = 0;
					ptrSimInfo->reportError("SamplingModule::sample_volume - Error: sampling did not converge in given number of iterations.");
				} else {
					ptrData->sv = i;
				}
			}

			return success;
		}
	}	// namespace sampling
}	// namespace pssalib
