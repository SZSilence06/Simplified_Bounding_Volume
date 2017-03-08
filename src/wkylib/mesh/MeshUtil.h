#ifndef WKY_MESH_UTIL_H
#define WKY_MESH_UTIL_H

#include "Common.h"

namespace WKYLIB
{
    namespace Mesh
    {
        void computeNormal(const matrixr_t& vertices, const matrixs_t& triangles, matrixr_t& normals);

        void computeCurveNormal2D(const matrixr_t& vertices, matrixr_t& normals);
    }
}

#endif
