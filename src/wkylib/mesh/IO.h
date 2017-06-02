#ifndef WKY_MESH_IO_H
#define WKY_MESH_IO_H

#include "Common.h"

namespace WKYLIB
{
    namespace Mesh
    {
        /**
         * @brief write points to a vtk format file.
         * @param file : output file path.
         * @param points : points to write. 
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool writePoints(const std::string &file, const std::vector<Eigen::Vector3d> &points);

        /**
         * @brief write 2d points to a vtk format file, with y coordinate set to 0, and z coordinate as y coordinate.
         * @param file : output file path.
         * @param points : points to write.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool writePoints2D(const std::string &file, const std::vector<Eigen::Vector2d> &points);

        /**
         * @brief read 2d curve from a obj format file, with y coordinate ignored, and z coordinate as y coordinate.
         * @param file : input file path.
         * @param vertices : matrix to store the vertex data;
         * @param lines : matrix to store the line data.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool readCurve2D(const std::string& file, std::vector<Eigen::Vector2d>& vertices,
                                 std::vector<Eigen::Vector2i>& lines);

        /**
         * @brief write 2D curve to a obj format file, with y coordinate set to 0, and z coordinate as y coordinate.
         * @param file : output file path.
         * @param vertices : vertices to output.
         * @param lines : lines to output.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool writeCurve2D(const std::string& file, const std::vector<Eigen::Vector2d>& vertices,
                                  const std::vector<Eigen::Vector2i>& lines);

        /**
         * @brief write triangle mesh to a obj format file.
         * @param file : output file path.
         * @param vertices : vertices to output.
         * @param triangles : triangles to output.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool writeMesh(const std::string &file, const std::vector<Eigen::Vector3d> &vertices,
                                 const std::vector<Eigen::Vector3i> &triangles);

        /**
         * @brief write 2D triangle mesh to a obj format file, with y coordinate set to 0, and z coordinate as y coordinate.
         * @param file : output file path.
         * @param vertices : vertices to output.
         * @param triangles : triangles to output.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool writeMesh2D(const std::string &file, const std::vector<Eigen::Vector2d> &vertices,
                                 const std::vector<Eigen::Vector3i> &triangles);

        /**
         * @brief write tetrahedron mesh to a vtk format file.
         * @param file : output file path.
         * @param vertices : vertices to output.
         * @param tetras : tetrahedrons to output.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool writeTetra(const std::string &file, const std::vector<Eigen::Vector3d> &vertices,
                                 const std::vector<Eigen::Vector4i> &tetras);
    }
}

#endif
