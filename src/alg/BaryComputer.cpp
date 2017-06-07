#include "BaryComputer.h"
#include "Shell.h"
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <iostream>

using namespace zjucad::matrix;

namespace SBV
{
  BaryComputer::BaryComputer(const matrixr_t &triangle)
    {
        buildInvA(triangle);
    }
    void BaryComputer::buildInvA(const matrixr_t &triangle)
    {
        this->invA = ones<double>(3, 3);
        this->invA(colon(0, 1), colon()) = triangle;
        if(inv(invA)) {
            std::cerr << "warning: degenerated triangle." << std::endl;
        }
    }
}
