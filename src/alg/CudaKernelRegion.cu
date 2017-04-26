#include "CudaKernelRegion.h"

namespace SBV
{
    CudaKernelRegion::CudaKernelRegion(const KernelRegion& cpu_kernel, const CudaPointer<CudaShell>& shell,
                                       const CudaPointer<CudaTriangulatedShell>& triangulation)
        : mShell(shell),
          mTriangulation(triangulation)
    {
        castMatToCuda(cpu_kernel.A, this->A);
        castMatToCuda_size_t(cpu_kernel.mLines, this->mLines);

        this->mPointType.assign(cpu_kernel.mPointType);
        this->mInvalidRegionType.assign(cpu_kernel.mInvalidRegionType);
        this->mClockwise.assign(cpu_kernel.mClockwise);

        std::vector<size_t> inner_samples;
        for(size_t sample : cpu_kernel.mInnerSamples)
        {
            inner_samples.push_back(sample);
        }
        mInnerSamples.assign(inner_samples);

        std::vector<size_t> outer_samples;
        for(size_t sample : cpu_kernel.mOuterSamples)
        {
            outer_samples.push_back(sample);
        }
        mOuterSamples.assign(outer_samples);
    }

    void CudaKernelRegion::castMatToCuda(const matrixr_t &matrix, CudaPointer<Eigen::MatrixXd> &cudaMat)
    {
        Eigen::MatrixXd eigenMat(matrix.size(1), matrix.size(2));
        for(int i = 0; i < matrix.size(1); i++)
        {
            for(int j = 0; j < matrix.size(2); j++)
            {
                eigenMat(i, j) = matrix(i, j);
            }
        }
        cudaMat.assign(eigenMat);
    }

    void CudaKernelRegion::castMatToCuda_size_t(const matrixs_t &matrix, CudaPointer<Eigen::MatrixXi> &cudaMat)
    {
        Eigen::MatrixXi eigenMat(matrix.size(1), matrix.size(2));
        for(int i = 0; i < matrix.size(1); i++)
        {
            for(int j = 0; j < matrix.size(2); j++)
            {
                eigenMat(i, j) = matrix(i, j);
            }
        }
        cudaMat.assign(eigenMat);
    }
}
