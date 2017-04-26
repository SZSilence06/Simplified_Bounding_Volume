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
        using CudaPointer = WKYLIB::Cuda::CudaPointer<T>;

        template<class T>
        using CudaVector = WKYLIB::Cuda::CudaVector<T>;

    public:
        CudaPointer<Eigen::MatrixXd> vertices;
        CudaPointer<Eigen::MatrixXd> triangles;
        CudaVector<PointType> vertType;
    };

    class CudaShell
    {
        template<class T>
        using CudaPointer = WKYLIB::Cuda::CudaPointer<T>;

        template<class T>
        using CudaVector = WKYLIB::Cuda::CudaVector<T>;

    public:
        CudaPointer<Eigen::MatrixXd> innerShell;
        CudaPointer<Eigen::MatrixXd> outerShell;
    };

    class CudaController
    {
    public:
        CudaController(const matrixr_t& innerShell, const matrixr_t& outerShell, const TriangulatedShell& triangulation);
        ~CudaController() = default;

    private:
        template<class T>
        using CudaPointer = WKYLIB::Cuda::CudaPointer<T>;

        void castMatToCuda(const matrixr_t& matrix, CudaPointer<Eigen::MatrixXd>& cudaMat);
        void castTriangulation(const TriangulatedShell& triangulation, CudaPointer<CudaTriangulatedShell>& cuda_triangulation);

    private:
        CudaPointer<CudaShell> mShell;
        CudaPointer<CudaTriangulatedShell> mCudaTriangulation;
    };
}

#endif
