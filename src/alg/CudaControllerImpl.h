#ifndef WKY_CUDA_CONTROLLER_IMPL_H
#define WKY_CUDA_CONTROLLER_IMPL_H

#include "Common.h"
#include "CudaShell.h"
#include "CudaTriangulatedShell.h"
#include "CudaKernelRegion.h"
#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>
#include "eigen3.3/Eigen/Dense"

namespace SBV
{
    class TriangulatedShell;
    class KernelRegion;

    class CudaControllerImpl
    {
    public:
        CudaControllerImpl();
        ~CudaControllerImpl() = default;

        void build(const matrixr_t& innerShell, const matrixr_t& outerShell, const TriangulatedShell& triangulation);

        void sample(double xmin, double xmax, double ymin, double ymax, double sampleRadius,
                    std::vector<matrixr_t>& output_samples);

        void buildKernelRegion(const KernelRegion& kernel);

        bool findCollapsePos_Boundary(bool isInner, const Eigen::Matrix3d& Q1, const Eigen::Matrix3d& Q2,
                                      matrixr_t& position, double& out_error);
        bool findCollapsePos_ZeroSet(const Eigen::Matrix3d& Q1, const Eigen::Matrix3d& Q2,
                                     matrixr_t& position, double& out_error);

    private:
        template<class T>
        using CudaPointer = WKYLIB::Cuda::CudaPointer<T>;

        template<class T>
        using CudaVector = WKYLIB::Cuda::CudaVector<T>;

        void castTriangulation(const TriangulatedShell& triangulation, CudaPointer<CudaTriangulatedShell>& cuda_triangulation);

    private:
        int mDeviceCount = 0;
        cudaDeviceProp mDeviceProperty;

        CudaPointer<CudaShell> mShell;
        CudaPointer<CudaTriangulatedShell> mTriangulation;
        CudaPointer<CudaKernelRegion> mKernel;
    };
}

#endif
