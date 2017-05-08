#ifndef WKY_SBV_COMMON_H
#define WKY_SBV_COMMON_H

#include <zjucad/matrix/matrix.h>

#pragma NVCC diagnostic ignored "-Wall"

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int
#if (defined __GNUC__) && (__GNUC__>4 || __GNUC_MINOR__>=7)
  #undef _GLIBCXX_ATOMIC_BUILTINS
  #undef _GLIBCXX_USE_INT128
#endif

namespace SBV
{
    using matrixr_t = zjucad::matrix::matrix<double>;
    using matrixs_t = zjucad::matrix::matrix<size_t>;

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
