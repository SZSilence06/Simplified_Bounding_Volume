#ifndef WKY_CUDA_TRIANGULATED_SHELL_H
#define WKY_CUDA_TRIANGULATED_SHELL_H

#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>
#include "eigen3.3/Eigen/Dense"
#include "PointType.h"

namespace SBV
{
    class CudaTriangulatedShell
    {
        template<class T>
        using CudaVector = WKYLIB::Cuda::CudaVector<T>;

    public:
        Eigen::MatrixXd vertices;
        Eigen::MatrixXd triangles;
        CudaVector<PointType> vertType;

    public:
        __device__ double getFValue(PointType pointType) const;
        __device__ double getFValue(size_t vert) const;
        __device__ double getSign(size_t vert) const;
    };
}

#endif
