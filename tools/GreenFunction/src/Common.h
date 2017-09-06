#ifndef GF_COMMON_H
#define GF_COMMON_H

#include <zjucad/matrix/matrix.h>
#include "matrix_mn.h"

using matrixr_t = zjucad::matrix::matrix<double>;
using matrixs_t = zjucad::matrix::matrix<size_t>;

typedef zjucad::matrix::matrix_mn<double, 2, 1> vec2_t;
typedef zjucad::matrix::matrix_mn<double, 3, 1> vec3_t;

#endif
