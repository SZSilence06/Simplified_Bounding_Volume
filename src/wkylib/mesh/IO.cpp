#include "IO.h"
#include <fstream>

namespace WKYLIB
{
    namespace Mesh
    {
        bool writePoints(const matrixr_t &points, const std::string &file)
        {
            std::ofstream out;

            out.open(file);
            if(out.fail())
            {
                return false;
            }

            out << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
                << points.size(2) << " float" << std::endl;
            for(int i = 0; i < points.size(2); i++)
            {
                out << points(0, i) << " " << points(1, i) << " " << points(2, i) << std::endl;
            }

            out << "CELLS " << points.size(2) << " " << points.size() * 2 << std::endl;
            for(int i = 0; i < points.size(2); i++)
            {
                out << "1 " << i << std::endl;
            }

            out << "CELL_TYPES " << points.size(2) << std::endl;
            for(int i = 0; i < points.size(2); i++)
            {
                out << "1" << std::endl;
            }

            out.close();

            return true;
        }
    }
}
