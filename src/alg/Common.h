#ifndef WKY_SBV_COMMON_H
#define WKY_SBV_COMMON_H

#include <zjucad/matrix/matrix.h>
#include "matrix_mn.h"

namespace SBV
{
    using matrixr_t = zjucad::matrix::matrix<double>;
    using matrixs_t = zjucad::matrix::matrix<size_t>;
  typedef zjucad::matrix::matrix_mn<double, 2, 1> vec2_t;
  typedef zjucad::matrix::matrix_mn<double, 3, 1> vec3_t;
  typedef zjucad::matrix::matrix_mn<double, 4, 1> vec4_t;
  typedef zjucad::matrix::matrix_mn<double, 2, 2> mat2x2_t;
  typedef zjucad::matrix::matrix_mn<double, 3, 3> mat3x3_t;
  typedef zjucad::matrix::matrix_mn<double, 4, 4> mat4x4_t;
    struct Mesh{
        matrixr_t vertices;
        matrixs_t triangles;
    };

    struct Curve{
        matrixr_t vertices;
        matrixs_t lines;
    };
}

#endif
