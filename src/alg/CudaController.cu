#include "CudaController.h"
#include <wkylib/Geometry/Util.h>

using namespace WKYLIB::Geometry;

namespace SBV
{
    __device__ double CudaTriangulatedShell::getFValue(PointType pointType) const
    {
        switch(pointType)
        {
        case POINT_BOUNDING_BOX:
            return 1;
        case POINT_OUTER:
            return 1;
        case POINT_INNER:
            return -1;
        }
        return 0;
    }

    __device__ double CudaTriangulatedShell::getFValue(size_t vert) const
    {
        return getFValue(vertType[vert]);
    }

    __device__ double CudaTriangulatedShell::getSign(size_t vert) const
    {
        double value = getFValue(vert);
        if(value > 0)
        {
            return 1;
        }
        else if(value == 0)
        {
            return 0;
        }
        return -1;
    }


    CudaController::CudaController(const matrixr_t &innerShell,
                                   const matrixr_t &outerShell,
                                   const TriangulatedShell &triangulation)
    {
        CudaShell shell;
        zju_mat_to_eigen(innerShell, shell.innerShell);
        zju_mat_to_eigen(outerShell, shell.outerShell);
        mShell.assign(shell);
        castTriangulation(triangulation, mCudaTriangulation);
    }

    void CudaController::castTriangulation(const TriangulatedShell &triangulation, CudaPointer<CudaTriangulatedShell> &cuda_triangulation)
    {
        cuda_triangulation.assign(CudaTriangulatedShell());
        zju_mat_to_eigen(triangulation.vertices, cuda_triangulation->vertices);
        zju_mat_to_eigen(triangulation.triangles, cuda_triangulation->triangles);
        cuda_triangulation->vertType.assign(triangulation.vertType);
    }
}
