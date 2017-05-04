#ifndef WKY_CUDA_KERNEL_REGION_H
#define WKY_CUDA_KERNEL_REGION_H

#include "Common.h"
#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>
#include "eigen3.3/Eigen/Dense"
#include "InvalidRegionType.h"
#include "CudaShell.h"
#include "CudaTriangulatedShell.h"

namespace SBV {
    class KernelRegion;

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


    private:
        void castMatToCuda(const matrixr_t& matrix, CudaPointer<Eigen::MatrixXd>& cudaMat);
        void castMatToCuda_size_t(const matrixs_t& matrix, CudaPointer<Eigen::MatrixXi>& cudaMat);

        __device__ bool isInvalidRegion(const Eigen::Vector2d& point) const;
        __device__ bool barycentric(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                                       const Eigen::Vector3d& p, Eigen::Vector3d& bary) const;
        __device__ bool barycentric_2D(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c,
                                       const Eigen::Vector2d& p, Eigen::Vector3d& bary) const;

    private:
        //these data are on gpu
        CudaPointer<Eigen::MatrixXi> mLines;
        CudaPointer<CudaShell> mShell;
        CudaPointer<CudaVector<size_t>> mInnerSamples;
        CudaPointer<CudaVector<size_t>> mOuterSamples;
        CudaPointer<CudaTriangulatedShell> mTriangulation;
        CudaPointer<PointType> mPointType;
        CudaPointer<InvalidRegionType> mInvalidRegionType;
        CudaPointer<bool> mClockwise;

        CudaPointer<Eigen::MatrixXd> A;
    };
}

#endif
