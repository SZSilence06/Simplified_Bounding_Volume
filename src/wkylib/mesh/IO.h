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
        WKY_API bool writePoints(const std::string& file, const matrixr_t& points);

        /**
         * @brief write 2d points to a vtk format file, with y coordinate set to 0, and z coordinate as y coordinate.
         * @param file : output file path.
         * @param points : points to write.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool writePoints2D(const std::string& file, const matrixr_t& points);

        /**
         * @brief read 2d curve from a obj format file, with y coordinate ignored, and z coordinate as y coordinate.
         * @param file : input file path.
         * @param vertices : matrix to store the vertex data;
         * @param lines : matrix to store the line data.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool readCurve2D(const std::string& file, matrixr_t& vertices, matrixs_t& lines);

        /**
         * @brief write 2D curve to a obj format file, with y coordinate set to 0, and z coordinate as y coordinate.
         * @param file : output file path.
         * @param vertices : vertices to output.
         * @param lines : lines to output.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool writeCurve2D(const std::string& file, const matrixr_t& vertices, const matrixs_t& lines);

        /**
         * @brief write 2D triangle mesh to a obj format file, with y coordinate set to 0, and z coordinate as y coordinate.
         * @param file : output file path.
         * @param vertices : vertices to output.
         * @param triangles : triangles to output.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool writeMesh2D(const std::string& file, const matrixr_t& vertices, const matrixs_t& triangles);

        /**
         * @brief write tetrahedral mesh to a vtk format file.
         * @param file : output file path.
         * @param vertices : vertices to output.
         * @param cells : tetrahedrons to output.
         * @return true if succeeded, and false otherwise.
         */
        WKY_API bool writeTetra(const std::string& file, const matrixr_t& vertices, const matrixs_t& cells);
    }
}

#endif
