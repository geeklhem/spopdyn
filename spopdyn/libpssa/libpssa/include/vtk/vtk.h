/**
 * \file vtk.h
 * \brief Declares classes for vtk visualization
 *
 * \details \c vtk class create vtk files format and call paraview,
 * for real time visualization of the species evolution
 *
 * \date   19.12.2012
 * \author Pietro Incardona
 */

#ifndef VTK_H_
#define VTK_H_

#include "../typedefs.h"

namespace pssalib
{
	namespace datamodel
	{
		class SimulationInfo;
	}
}

namespace pssalib
{
	namespace vtk
	{
		////////////////////////////////////////
		//! Write a sequence vtk file
		class vtk
		{
			private:

				std::string file;
				int sequence;
			  
			public:

				//! Create a sequence of VTK file
				void WriteVtkFile(std::string str, pssalib::datamodel::SimulationInfo * dt, int sq);
				void Visualize();
				
				vtk()	{};
				~vtk()	{};
		};
	}
}

#endif /* VTK_H_ */
