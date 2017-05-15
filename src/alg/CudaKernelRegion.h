#ifndef WKY_CUDA_KERNEL_REGION_H
#define WKY_CUDA_KERNEL_REGION_H

#include "Common.h"
#include <wkylib/Cuda/CudaEigen.h>
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

        __host__ __device__ void printTest();
        
        inline __host__ __device__ const CudaVector<size_t>& getInnerSamples() const { return this->mInnerSamples; }
        inline __host__ __device__ const CudaVector<size_t>& getOuterSamples() const { return this->mOuterSamples; }
        inline __host__ __device__ const CudaPointer<CudaShell>& getShell() const { return this->mShell; }

    private:
        __device__ bool isInvalidRegion(const Eigen::Vector2d& point) const;
        __device__ bool barycentric(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                                       const Eigen::Vector3d& p, Eigen::Vector3d& bary) const;
        __device__ bool barycentric_2D(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c,
                                       const Eigen::Vector2d& p, Eigen::Vector3d& bary) const;

        __host__ __device__ void printEigenMat(const CudaPointer<Eigen::MatrixXd>& mat, const char* name);
        __host__ __device__ void printEigenMat(const CudaPointer<Eigen::MatrixXi>& mat, const char* name);
        __host__ __device__ void printVector(const CudaVector<size_t>& vec, const char* name);
        __host__ __device__ void printVector(const CudaVector<PointType>& vec, const char* name);

    private:
        //these data are on gpu
        CudaVector<Eigen::Vector2i> mLines;
        CudaPointer<CudaShell> mShell;
        CudaPointer<CudaTriangulatedShell> mTriangulation;
        CudaVector<size_t> mInnerSamples;
        CudaVector<size_t> mOuterSamples;
        PointType mPointType;
        InvalidRegionType mInvalidRegionType;
        bool mClockwise;

        CudaPointer<Eigen::MatrixXd> A;

        friend class CudaControllerImpl;
    };
}

#endif
