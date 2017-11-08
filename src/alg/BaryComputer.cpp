#include "BaryComputer.h"
#include "Shell.h"
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <iostream>

using namespace zjucad::matrix;

namespace SBV
{
    BaryComputer::BaryComputer(const matrixr_t &tetra)
    {
        buildInvA(tetra);
    }
    void BaryComputer::buildInvA(const matrixr_t &tetra)
    {
        this->invA = ones<double>(4, 4);
        this->invA(colon(0, 2), colon()) = tetra;
        if(inv(invA)) {
            std::cerr << "warning: degenerated triangle. line " << __LINE__ << std::endl;
        }
    }
}
