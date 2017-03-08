#ifndef WKY_SBV_COMMON_H
#define WKY_SBV_COMMON_H

#include <zjucad/matrix/matrix.h>

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
