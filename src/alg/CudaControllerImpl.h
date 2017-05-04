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
        CudaControllerImpl() = default;
        ~CudaControllerImpl() = default;

        void build(const matrixr_t& innerShell, const matrixr_t& outerShell, const TriangulatedShell& triangulation);

        void sample(double xmin, double xmax, double ymin, double ymax, double sampleRadius,
                    std::vector<matrixr_t>& output_samples);

        void buildKernelRegion(const KernelRegion& kernel);

    private:
        template<class T>
        using CudaPointer = WKYLIB::Cuda::CudaPointer<T>;

        void castTriangulation(const TriangulatedShell& triangulation, CudaPointer<CudaTriangulatedShell>& cuda_triangulation);

    private:
        CudaPointer<CudaShell> mShell;
        CudaPointer<CudaTriangulatedShell> mTriangulation;
        CudaPointer<CudaKernelRegion> mKernel;
    };
}

#endif
