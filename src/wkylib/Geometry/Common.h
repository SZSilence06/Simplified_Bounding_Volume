#ifndef WKY_GEOMETRY_COMMON_H
#define WKY_GEOMETRY_COMMON_H

#include "matrix_mn.h"
#include <wkylib/Common.h>
#include <zjucad/matrix/matrix.h>

namespace WKYLIB
{
    namespace Geometry
    {
        using matrixr_t = zjucad::matrix::matrix<double>;
        using matrixs_t = zjucad::matrix::matrix<size_t>;
        typedef zjucad::matrix::matrix_mn<double, 2, 1> vec2_t;
        typedef zjucad::matrix::matrix_mn<double, 3, 1> vec3_t;
        typedef zjucad::matrix::matrix_mn<double, 4, 1> vec4_t;
        typedef zjucad::matrix::matrix_mn<double, 2, 2> mat2x2_t;
        typedef zjucad::matrix::matrix_mn<double, 3, 3> mat3x3_t;
        typedef zjucad::matrix::matrix_mn<double, 4, 4> mat4x4_t;
    }
}

#endif
