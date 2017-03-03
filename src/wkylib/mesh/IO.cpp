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

            for(int i = 0; i < points.size(2); i++)
            {
                out << "v " << points(0, i) << " " << points(1, i) << " " << points(2, i) << std::endl;
            }

            out.close();

            return true;
        }
    }
}
