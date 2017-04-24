#include "CudaController.h"

namespace SBV
{
    CudaController::CudaController(const matrixr_t &innerShell,
                                   const matrixr_t &outerShell,
                                   const TriangulatedShell &triangulation)
    {
        castMatToCuda(innerShell, mInnerShell);
        castMatToCuda(outerShell, mOuterShell);
    }

    void CudaController::castMatToCuda(const matrixr_t &matrix, CudaPointer<Eigen::MatrixXd> &cudaMat)
    {
        Eigen::MatrixXd eigenMat(matrix.size(1), matrix.size(2));
        for(int i = 0; i < matrix.size(1); i++)
        {
            for(int j = 0; j < matrix.size(2); j++)
            {
                eigenMat(i, j) = matrix(i, j);
            }
        }
        cudaMat = eigenMat;
    }
}
