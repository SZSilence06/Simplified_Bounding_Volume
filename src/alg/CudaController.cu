#include "CudaController.h"

namespace SBV
{
    CudaController::CudaController(const matrixr_t &innerShell,
                                   const matrixr_t &outerShell,
                                   const TriangulatedShell &triangulation)
    {
        CudaShell shell;
        castMatToCuda(innerShell, shell.innerShell);
        castMatToCuda(outerShell, shell.outerShell);
        mShell.assign(shell);
        castTriangulation(triangulation, mCudaTriangulation);
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
        cudaMat.assign(eigenMat);
    }

    void CudaController::castTriangulation(const TriangulatedShell &triangulation, CudaPointer<CudaTriangulatedShell> &cuda_triangulation)
    {
        castMatToCuda(triangulation.vertices, cuda_triangulation->vertices);
        castMatToCuda(triangulation.triangles, cuda_triangulation->triangles);
        cuda_triangulation->vertType.assign(triangulation.vertType);
    }
}
