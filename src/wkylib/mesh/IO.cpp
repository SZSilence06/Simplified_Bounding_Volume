#include "IO.h"
#include <fstream>
#include <sstream>

using namespace zjucad::matrix;

namespace WKYLIB
{
    namespace Mesh
    {
        bool writePoints(const std::string &file, const matrixr_t &points)
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

        bool writePoints2D(const std::string &file, const matrixr_t &points)
        {
            if(points.size(1) != 2)
            {
                throw std::invalid_argument("You can only use 2d matrices in writePoints2D().");
            }

            matrixr_t points_3D(3, points.size(2));
            for(int i = 0; i < points.size(2); i++)
            {
                points_3D(0, i) = points(0, i);
                points_3D(1, i) = 0;
                points_3D(2, i) = points(1, i);
            }
            return writePoints(file, points_3D);
        }


        bool readCurve2D(const std::string& file, matrixr_t& vertices, matrixs_t& lines)
        {
            std::ifstream in;

            in.open(file);
            if(in.fail())
            {
                return false;
            }

            std::vector<matrixr_t> curve;
            std::vector<matrixs_t> seg;
            while(in)
            {
                char buf[1024];
                in.getline(buf, sizeof(buf));

                std::stringstream ss;
                ss << buf;
                char sign;
                ss >> sign;
                if(sign == 'v')
                {
                    double a, b, c;
                    ss >> a >> b >> c;

                    matrixr_t vert(2,1);
                    vert[0] = a;
                    vert[1] = c;
                    curve.push_back(vert);
                }
                else if(sign == 'l')
                {
                    int a, b;
                    ss >> a >> b;

                    matrixs_t line(2, 1);
                    line[0] = a - 1;
                    line[1] = b - 1;
                    seg.push_back(line);
                }
            }

            vertices.resize(2, curve.size());
            for(int i = 0; i < curve.size(); i++)
            {
                vertices(colon(), i) = curve[i];
            }

            lines.resize(2, seg.size());
            for(int i = 0; i < seg.size(); i++)
            {
                lines(colon(), i) = seg[i];
            }

            return true;
        }

        bool writeCurve2D(const std::string& file, const matrixr_t& vertices, const matrixs_t& lines)
        {
            if(vertices.size(1) != 2 || lines.size(1) != 2)
            {
                throw std::invalid_argument("You can only use 2d matrices in writeCurve2D().");
            }

            std::ofstream out;

            out.open(file);
            if(out.fail())
            {
                return false;
            }

            for(int i = 0; i < vertices.size(2); i++)
            {
                out << "v " << vertices(0, i) << " 0 " << vertices(1, i) << std::endl;
            }
            for(int i = 0; i < lines.size(2); i++)
            {
                out << "l " << lines(0, i) + 1 << " " << lines(1, i) + 1 << std::endl;
            }

            out.close();

            return true;
        }

        bool writeMesh2D(const std::string &file, const matrixr_t &vertices, const matrixs_t &triangles)
        {
            if(vertices.size(1) != 2 || triangles.size(1) != 3)
            {
                throw std::invalid_argument("You can only use 2d matrix for 'vertices' and 3d matrix for 'trianlges' in writeMesh2D().");
            }

            std::ofstream out;

            out.open(file);
            if(out.fail())
            {
                return false;
            }

            for(int i = 0; i < vertices.size(2); i++)
            {
                out << "v " << vertices(0, i) << " 0 " << vertices(1, i) << std::endl;
            }
            for(int i = 0; i < triangles.size(2); i++)
            {
                out << "f " << triangles(0, i) + 1 << " " << triangles(1, i) + 1 << " " << triangles(2, i) + 1 << std::endl;
            }

            out.close();

            return true;
        }
    }
}
