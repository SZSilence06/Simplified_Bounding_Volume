#include "CudaTriangulatedShell.h"

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
}
