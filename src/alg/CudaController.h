#ifndef WKY_CUDA_CONTROLLER_H
#define WKY_CUDA_CONTROLLER_H

#include "TriangulatedShell.h"
#include <thrust/device_vector.h>
#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>
#include "eigen3.3/Eigen/Dense"

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

    class CudaShell
    {
    public:
        Eigen::MatrixXd innerShell;
        Eigen::MatrixXd outerShell;
    };

    class CudaController
    {
    public:
        CudaController(const matrixr_t& innerShell, const matrixr_t& outerShell, const TriangulatedShell& triangulation);
        ~CudaController() = default;

    private:
        template<class T>
        using CudaPointer = WKYLIB::Cuda::CudaPointer<T>;

        void castTriangulation(const TriangulatedShell& triangulation, CudaPointer<CudaTriangulatedShell>& cuda_triangulation);

    private:
        CudaPointer<CudaShell> mShell;
        CudaPointer<CudaTriangulatedShell> mCudaTriangulation;
    };
}

#endif
