#ifndef WKY_MESH_COMMON_H
#define WKY_MESH_COMMON_H

#include <zjucad/matrix/matrix.h>

namespace WKYLIB
{
    namespace Mesh
    {
        using matrixr_t = zjucad::matrix::matrix<double>;
        using matrixs_t = zjucad::matrix::matrix<size_t>;
    }
}

#endif
