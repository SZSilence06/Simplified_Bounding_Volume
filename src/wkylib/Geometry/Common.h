#ifndef WKY_GEOMETRY_COMMON_H
#define WKY_GEOMETRY_COMMON_H

#include <wkylib/Common.h>

#include <zjucad/matrix/matrix.h>

namespace WKYLIB
{
    namespace Geometry
    {
        using matrixr_t = zjucad::matrix::matrix<double>;
        using matrixs_t = zjucad::matrix::matrix<size_t>;
    }
}

#endif
