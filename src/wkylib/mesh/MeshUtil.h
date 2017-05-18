#ifndef WKY_MESH_UTIL_H
#define WKY_MESH_UTIL_H

#include "Common.h"

namespace WKYLIB
{
    namespace Mesh
    {
        WKY_API void computeNormal(const matrixr_t& vertices, const matrixs_t& triangles, matrixr_t& normals);

        WKY_API void computeNormal2D(const std::vector<Eigen::Vector2d>& vertices, const std::vector<Eigen::Vector2i>&lines,
                                           std::vector<Eigen::Vector2d>& normals);
    }
}

#endif
