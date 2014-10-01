#include "../../include/vtk/vtk.h"
#include <fstream>
#include "../../include/datamodel/SimulationInfo.h"
#include <sstream>

 namespace pssalib
{
	namespace vtk
	{
		void vtk::WriteVtkFile(std::string str, pssalib::datamodel::SimulationInfo * dt, int sq)
		{
			// Open a file

			std::ofstream f;
			
			std::stringstream sstm;
			sstm << str.c_str() << sq << ".vtk";
			str = sstm.str();
			
			f.open(str.c_str(),std::fstream::out);
			
			if (f.is_open() == false)
			{std::cerr << "Error: cannot open vtk output file\n"; return;}
			
			// VTK header file

			f << "# vtk DataFile Version 3.0\n";
			f << "Partial propensity vtk frame: " << sq << " volumes: " << dt->getDataModel()->unSubvolumes << "\n";

			f << "ASCII\n";
			f << "DATASET STRUCTURED_POINTS\n";
			if (dt->eVolumeType == pssalib::datamodel::SimulationInfo::VT_Homogeneous2d)
			{f << "DIMENSIONS " << dt->unNumGridPoints+1 << " " << dt->unNumGridPoints+1 << " " << 1 << "\n";}
			else
			{f << "DIMENSIONS " << dt->unNumGridPoints+1 << " " << dt->unNumGridPoints+1 << " " << dt->unNumGridPoints+1 << "\n";}
			f << "ORIGIN 0 0 0\n";
			f << "SPACING 1 1 1\n";
			f << "CELL_DATA " << dt->getDataModel()->unSubvolumes << "\n";
			f << "SCALARS population int\n";
			f << "LOOKUP_TABLE default\n";

			for (int i = 0 ; i < dt->getDataModel()->unSubvolumes ; i++)
			{
				f << dt->getDataModel()->arSubvolumes[i].aruN[1] << "\n";
			}
			
			f << "SCALARS population2 int\n";
			f << "LOOKUP_TABLE default\n";

			for (int i = 0 ; i < dt->getDataModel()->unSubvolumes ; i++)
			{
				f << dt->getDataModel()->arSubvolumes[i].aruN[2] << "\n";
			}
			
			f.close();
		}
		
		void vtk::Visualize()
		{
			
		}
	}
}

