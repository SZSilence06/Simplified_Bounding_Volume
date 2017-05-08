#ifndef WKY_CUDA_SAMPLING_TREE_H
#define WKY_CUDA_SAMPLING_TREE_H

#include "Common.h"
#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>
#include "eigen3.3/Eigen/Dense"

namespace SBV
{
    class CudaKernelRegion;

    class CudaSamplingTree
    {
    private:
         template<class T>
         using CudaPointer = WKYLIB::Cuda::CudaPointer<T>;

         template<class T>
         using CudaVector = WKYLIB::Cuda::CudaVector<T>;

    public:
        CudaSamplingTree(CudaPointer<CudaKernelRegion> kernel, double xmax, double xmin, double ymax, double ymin, double sampleRadius);

        inline const CudaVector<Eigen::Vector2d>& getSamples() const { return *mSamples; }

    private:
        void sample(double xmin, double xmax, double ymin, double ymax);

    private:
        CudaPointer<CudaKernelRegion> mKernel;
        double mSampleRadius;
        CudaPointer<CudaVector<Eigen::Vector2d>> mSamples;
    };
}

#endif
