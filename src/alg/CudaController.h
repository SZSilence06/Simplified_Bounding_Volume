#ifndef WKY_CUDA_CONTROLLER_H
#define WKY_CUDA_CONTROLLER_H

#include "TriangulatedShell.h"
#include <thrust/device_vector.h>
#include <wkylib/Cuda/CudaPointer.h>
#include "eigen3.3/Eigen/Dense"

namespace SBV
{
    class CudaController
    {
    public:
        CudaController(const matrixr_t& innerShell, const matrixr_t& outerShell, const TriangulatedShell& triangulation);
        ~CudaController() = default;

    private:
        template<class T>
        using CudaPointer = typename WKYLIB::Cuda::CudaPointer<T>;

        void castMatToCuda(const matrixr_t& matrix, CudaPointer<Eigen::MatrixXd>& cudaMat);

    private:
        CudaPointer<Eigen::MatrixXd> mInnerShell;
        CudaPointer<Eigen::MatrixXd> mOuterShell;
        TriangulatedShell* mCudaTriangulation = nullptr;
    };
}

#endif
