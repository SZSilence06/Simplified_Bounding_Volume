#include "IO.h"
#include <fstream>
#include <sstream>

using namespace zjucad::matrix;

namespace WKYLIB
{
    namespace Mesh
    {
        bool writePoints(const std::string &file, const std::vector<Eigen::Vector3d> &points)
        {
            std::ofstream out;

            out.open(file);
            if(out.fail())
            {
                return false;
            }

            out << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
                << points.size() << " float" << std::endl;
            for(size_t i = 0; i < points.size(); i++)
            {
                out << points[i][0] << " " << points[i][1] << " " << points[i][2] << std::endl;
            }

            out << "CELLS " << points.size() << " " << points.size() * 2 << std::endl;
            for(size_t i = 0; i < points.size(); i++)
            {
                out << "1 " << i << std::endl;
            }

            out << "CELL_TYPES " << points.size() << std::endl;
            for(size_t i = 0; i < points.size(); i++)
            {
                out << "1" << std::endl;
            }

            out.close();

            return true;
        }

        bool writePoints2D(const std::string &file, const std::vector<Eigen::Vector2d> &points)
        {
            std::vector<Eigen::Vector3d> points_3D(points.size());
            for(size_t i = 0; i < points_3D.size(); i++)
            {
                points_3D[i][0] = points[i][0];
                points_3D[i][1] = 0;
                points_3D[i][2] = points[i][1];
            }
            return writePoints(file, points_3D);
        }


        bool readCurve2D(const std::string& file, std::vector<Eigen::Vector2d>& vertices,
                         std::vector<Eigen::Vector2i>& lines)
        {
            std::ifstream in;

            in.open(file);
            if(in.fail())
            {
                return false;
            }

            vertices.clear();
            lines.clear();
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

                    Eigen::Vector2d vert;
                    vert[0] = a;
                    vert[1] = c;
                    vertices.push_back(vert);
                }
                else if(sign == 'l')
                {
                    int a, b;
                    ss >> a >> b;

                    Eigen::Vector2i line;
                    line[0] = a - 1;
                    line[1] = b - 1;
                    lines.push_back(line);
                }
            }
            return true;
        }

        bool writeCurve2D(const std::string& file, const std::vector<Eigen::Vector2d>& vertices,
                          const std::vector<Eigen::Vector2i>& lines)
        {
            std::ofstream out;

            out.open(file);
            if(out.fail())
            {
                return false;
            }

            for(size_t i = 0; i < vertices.size(); i++)
            {
                out << "v " << vertices[i][0] << " 0 " << vertices[i][1] << std::endl;
            }
            for(size_t i = 0; i < lines.size(); i++)
            {
                out << "l " << lines[i][0] + 1 << " " << lines[i][1] + 1 << std::endl;
            }

            out.close();

            return true;
        }

        bool writeMesh(const std::string &file, const std::vector<Eigen::Vector3d> &vertices,
                       const std::vector<Eigen::Vector3i> &triangles)
        {
            std::ofstream out;

            out.open(file);
            if(out.fail())
            {
                return false;
            }

            for(size_t i = 0; i < vertices.size(); i++)
            {
                out << "v " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << std::endl;
            }
            for(size_t i = 0; i < triangles.size(); i++)
            {
                out << "f " << triangles[i][0] + 1 << " " << triangles[i][1] + 1 << " " << triangles[i][2] + 1 << std::endl;
            }

            out.close();

            return true;
        }

        bool writeMesh2D(const std::string &file, const std::vector<Eigen::Vector2d> &vertices,
                         const std::vector<Eigen::Vector3i> &triangles)
        {
            std::ofstream out;

            out.open(file);
            if(out.fail())
            {
                return false;
            }

            for(size_t i = 0; i < vertices.size(); i++)
            {
                out << "v " << vertices[i][0] << " 0 " << vertices[i][1] << std::endl;
            }
            for(size_t i = 0; i < triangles.size(); i++)
            {
                out << "f " << triangles[i][0] + 1 << " " << triangles[i][1] + 1 << " " << triangles[i][2] + 1 << std::endl;
            }

            out.close();

            return true;
        }

        bool writeTetra(const std::string &file, const std::vector<Eigen::Vector3d> &vertices,
                        const std::vector<Eigen::Vector4i> &tetras)
        {
            std::ofstream out;

            out.open(file);
            if(out.fail())
            {
                return false;
            }

            out << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
                << vertices.size() << " float" << std::endl;
            for(size_t i = 0; i < vertices.size(); i++)
            {
                out << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << std::endl;
            }

            out << "CELLS " << tetras.size() << " " << tetras.size() * 5 << std::endl;
            for(size_t i = 0; i < tetras.size(); i++)
            {
                out << "4";
                for(int j = 0; j < 4; j++)
                {
                    out << " "<< tetras[i][j];
                }
                out << std::endl;
            }

            out << "CELL_TYPES " << tetras.size() << std::endl;
            for(size_t i = 0; i < tetras.size(); i++)
            {
                out << "10" << std::endl;
            }

            out.close();
            return true;
        }
    }
}
