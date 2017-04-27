#ifndef WKY_CUDA_KERNEL_REGION_H
#define WKY_CUDA_KERNEL_REGION_H

#include "TriangulatedShell.h"
#include <thrust/device_vector.h>
#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>
#include "eigen3.3/Eigen/Dense"
#include "KernelRegion.h"
#include "CudaController.h"

namespace SBV {
    class CudaKernelRegion
    {
        template<class T>
        using CudaPointer = WKYLIB::Cuda::CudaPointer<T>;

        template<class T>
        using CudaVector = WKYLIB::Cuda::CudaVector<T>;

    public:
        CudaKernelRegion(const KernelRegion& cpu_kernel, const CudaPointer<CudaShell>& shell,
                         const CudaPointer<CudaTriangulatedShell>& triangulation);

        __device__ bool contains(const Eigen::Vector2d& point) const;
        __device__ bool IsInvalidRegion(const Eigen::Vector2d& point) const;

    private:
        void castMatToCuda(const matrixr_t& matrix, CudaPointer<Eigen::MatrixXd>& cudaMat);
        void castMatToCuda_size_t(const matrixs_t& matrix, CudaPointer<Eigen::MatrixXi>& cudaMat);

    private:
        //these data are on gpu
        CudaPointer<Eigen::MatrixXi> mLines;
        CudaPointer<CudaShell> mShell;
        CudaVector<size_t> mInnerSamples;
        CudaVector<size_t> mOuterSamples;
        CudaPointer<CudaTriangulatedShell> mTriangulation;
        CudaPointer<PointType> mPointType;
        CudaPointer<InvalidRegionType> mInvalidRegionType;
        CudaPointer<bool> mClockwise;

        CudaPointer<Eigen::MatrixXd> A;
    };
}

#endif
