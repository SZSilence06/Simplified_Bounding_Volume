#ifndef WKY_MESH_UTIL_H
#define WKY_MESH_UTIL_H

#include "Common.h"

namespace WKYLIB
{
    namespace Mesh
    {
        WKY_API void computeNormal(const matrixr_t& vertices, const matrixs_t& triangles, matrixr_t& normals);

        WKY_API void computeNormal2D(const matrixr_t& vertices, const matrixr_t& lines, matrixr_t& normals);
    }
}

#endif
