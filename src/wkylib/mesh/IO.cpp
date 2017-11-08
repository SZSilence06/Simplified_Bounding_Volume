#include "IO.h"
#include <fstream>
#include <sstream>

using namespace zjucad::matrix;

namespace WKYLIB
{
    namespace Mesh
    {
        static bool writePoints_vtk(const std::string &file, const matrixr_t &points)
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

            out << "CELLS " << points.size(2) << " " << points.size(2) * 2 << std::endl;
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

        static bool writePoints_pcd(const std::string &file, const matrixr_t &points)
        {
            std::ofstream out;
            out.open(file);
            if(out.fail())
            {
                return false;
            }

            out << "VERSION .7" << std::endl;
            out << "FIELDS x y z" << std::endl;
            out << "SIZE 4 4 4" << std::endl;
            out << "TYPE F F F" << std::endl;
            out << "COUNT 1 1 1" << std::endl;
            out << "WIDTH " << points.size(2) << std::endl;
            out << "HEIGHT 1" << std::endl;
            out << "VIEWPOINT 0 0 0 1 0 0 0" << std::endl;
            out << "POINTS " << points.size(2) << std::endl;
            out << "DATA ascii" << std::endl;

            for(int i = 0; i < points.size(2); i++)
            {
                out << points(0, i) << " " << points(1, i) << " " << points(2, i) << std::endl;
            }

            out.close();

            return true;
        }

        static bool writePoints_obj(const std::string &file, const matrixr_t &points)
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

        bool writePoints(const std::string &file, const matrixr_t &points)
        {
            std::string ext = file.substr(file.find_last_of('.') + 1);
            if(ext == "vtk")
            {
                return writePoints_vtk(file, points);
            }
            if(ext == "pcd")
            {
                return writePoints_pcd(file, points);
            }
            if(ext == "obj")
            {
                return writePoints_obj(file, points);
            }
            return false;
        }

        static bool readPoints_vtk(const std::string& file, matrixr_t& points)
        {
            std::ifstream in;


            in.open(file);
            if(in.fail())
            {
                return false;
            }

            std::string head;

            while(in)
            {
                char buf[1024];
                in.getline(buf, sizeof(buf));
                std::stringstream ss;
                ss << buf;

                ss >> head;
                if(head == "POINTS")
                {
                    int point_number;
                    ss >> point_number;
                    points.resize(3, point_number);
                    for(int i = 0; i < point_number; i++)
                    {
                        double x, y, z;
                        std::stringstream ss2;
                        in.getline(buf, sizeof(buf));
                        ss2 << buf;
                        ss2 >> x >> y >> z;
                        points(0, i) = x;
                        points(1, i) = y;
                        points(2, i) = z;
                    }
                    return true;
                }
            }

            return false;
        }

        bool readPoints(const std::string &file, matrixr_t &points)
        {
            std::string ext = file.substr(file.find_last_of('.') + 1);
            if(ext == "vtk")
            {
                return readPoints_vtk(file, points);
            }
            return false;
        }

        static bool writePointsAndNormals_obj(const std::string &file, const matrixr_t &points, const matrixr_t& normals)
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
            for(int i = 0; i < normals.size(2); i++)
            {
                out << "vn " << normals(0, i) << " " << normals(1, i) << " " << normals(2, i) << std::endl;
            }

            out.close();

            return true;
        }

        static bool writePointsAndNormals_vtk(const std::string &file, const matrixr_t &points, const matrixr_t& normals)
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

            out << "CELLS " << points.size(2) << " " << points.size(2) * 2 << std::endl;
            for(int i = 0; i < points.size(2); i++)
            {
                out << "1 " << i << std::endl;
            }

            out << "CELL_TYPES " << points.size(2) << std::endl;
            for(int i = 0; i < points.size(2); i++)
            {
                out << "1" << std::endl;
            }

            out << "CELL_DATA " << points.size(2) << std::endl;
            out << "NORMALS normals float" << std::endl;
            for(int i = 0; i < normals.size(2); i++)
            {
                out << normals(0, i) << " " << normals(1, i) << " " << normals(2, i) << std::endl;
            }

            out.close();

            return true;
        }

        bool writePointsAndNormals(const std::string &file, const matrixr_t &points, const matrixr_t& normals)
        {
            std::string ext = file.substr(file.find_last_of('.') + 1);
            if(ext == "obj")
            {
                return writePointsAndNormals_obj(file, points, normals);
            }
            else if(ext == "vtk")
            {
                return writePointsAndNormals_vtk(file, points, normals);
            }
            return false;
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
            char sign = 0;

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
                sign = 0;
                std::stringstream ss;
                ss << buf;

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

        bool writeTetra(const std::string& file, const matrixr_t& vertices, const matrixs_t& cells)
        {
            std::ofstream out;

            out.open(file);
            if(out.fail())
            {
                return false;
            }

            out << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
                << vertices.size(2) << " float" << std::endl;
            for(int i = 0; i < vertices.size(2); i++)
            {
                out << vertices(0, i) << " " << vertices(1, i) << " " << vertices(2, i) << std::endl;
            }

            out << "CELLS " << cells.size(2) << " " << cells.size(2) * 5 << std::endl;
            for(int i = 0; i < cells.size(2); i++)
            {
                out << "4";
                for(int j = 0; j < 4; j++)
                {
                    out << " " << cells(j, i);
                }
                out << std::endl;
            }

            out << "CELL_TYPES " << cells.size(2) << std::endl;
            for(int i = 0; i < cells.size(2); i++)
            {
                out << "10" << std::endl;
            }

            out.close();

            return true;
        }

        WKY_API bool writeMeshAndNormals(const std::string& file, const matrixr_t& vertices, const matrixs_t& triangles, const matrixr_t& normals)
        {
            if(vertices.size(1) != 3 || triangles.size(1) != 3 || normals.size(1) != 3)
            {
                throw std::invalid_argument("You can only use 3d matrix for 'vertices' and 3d matrix for 'trianlges' and 3d matrix for 'normals' in writeMeshAndNormals().");
            }

            std::ofstream out;
            out.open(file);
            if(out.fail())
            {
                return false;
            }

            for(int i = 0; i < vertices.size(2); i++)
            {
                out << "v " << vertices(0, i) << " " << vertices(1, i) << " " << vertices(2, i) << std::endl;
            }
            for(int i = 0; i < normals.size(2); i++)
            {
                out << "vn " << normals(0, i) << " " << normals(1, i) << " " << normals(2, i) << std::endl;
            }
            for(int i = 0; i < triangles.size(2); i++)
            {
                out << "f " << triangles(0, i) + 1 << "//" << triangles(0, i) + 1 << " "
                    << triangles(1, i) + 1 << "//" << triangles(1, i) + 1 << " "
                    << triangles(2, i) + 1 << "//" << triangles(2, i) + 1 << std::endl;
            }

            out.close();
            return true;
        }
    }
}
