#ifndef WKY_MESH_UTIL_H
#define WKY_MESH_UTIL_H

#include "Common.h"

namespace WKYLIB
{
    namespace Mesh
    {
        WKY_API void computeNormal(const std::vector<Eigen::Vector3d>& vertices, const std::vector<Eigen::Vector3i>&triangles,
                                   std::vector<Eigen::Vector3d>& normals);

        WKY_API void computeNormal2D(const std::vector<Eigen::Vector2d>& vertices, const std::vector<Eigen::Vector2i>&lines,
                                           std::vector<Eigen::Vector2d>& normals);
    }
}

#endif
