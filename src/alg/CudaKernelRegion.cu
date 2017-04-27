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

    /*__device__ bool CudaKernelRegion::contains(const Eigen::Vector2d &point) const
    {
        Eigen::Vector3d homo;
        homo[0] = point[0];
        homo[1] = point[1];
        homo[2] = 1;
        Eigen::VectorXd result = A * homo;
        for(int i = 0; i < result.rows(); i++)
        {
            if(result[i] > 0)
            {
                return false;
            }
        }
        if(isInvalidRegion(point))
        {
            return false;
        }
        return true;
    }

    __device__ bool CudaKernelRegion::IsInvalidRegion(const Eigen::Vector2d& point) const
    {
        auto& lines = *mLines;
        auto& vertices = *mTriangulation->vertices;
        for(int i = 0; i < lines.rows(); i++)
        {
            Eigen::MatrixXd triangle(2, 3);
            triangle.col(0) = vertices.col(lines(0, i));
            triangle.col(1) = vertices.col(lines(1, i));
            triangle.col(2) = point;

            for(int i = 0; i < mInnerSamples.size(); i++)
            {
                size_t sample = mInnerSamples.getElements()[i];
                matrixr_t bary;
                if(WKYLIB::barycentric_2D(mShell.mInnerShell(colon(), sample), triangle, bary))
                {
                    //the point is inside the tetrahedron
                    double f0 = mTriangulation.getFValue(mLines(0, i));
                    double f1 = mTriangulation.getFValue(mLines(1, i));
                    double f2 = mTriangulation.getFValue(mPointType);

                    double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2];
                    if(f > 0)
                    {
                        return true;
                    }
                }
            }

            for(size_t sample : mOuterSamples)
            {
                matrixr_t bary;
                if(WKYLIB::barycentric_2D(mShell.mOuterShell(colon(), sample), triangle, bary))
                {
                    //the point is inside the tetrahedron
                    double f0 = mTriangulation.getFValue(mLines(0, i));
                    double f1 = mTriangulation.getFValue(mLines(1, i));
                    double f2 = mTriangulation.getFValue(mPointType);

                    double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2];
                    if(f < 0)
                    {
                        return true;
                    }
                }
            }
        }

        return false;
    }*/
}
