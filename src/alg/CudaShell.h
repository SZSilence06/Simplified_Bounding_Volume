#ifndef WKY_CUDA_SHELL_H
#define WKY_CUDA_SHELL_H

#include "Common.h"
#include "eigen3/Eigen/Dense"
#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>
#include <wkylib/Cuda/CudaEigen.h>

namespace SBV
{
    class CudaShell
    {
        template<class T>
        using CudaPointer = WKYLIB::Cuda::CudaPointer<T>;

        template<class T>
        using CudaVector = WKYLIB::Cuda::CudaVector<T>;

    public:
        CudaVector<Point> innerShell;
        CudaVector<Point> outerShell;
    };
}

#endif
