#include "CudaKernelRegion.h"
#include "Common.h"
#include "KernelRegion.h"
#include <wkylib/Geometry/Util.h>

using namespace WKYLIB::Geometry;

namespace SBV
{
    CudaKernelRegion::CudaKernelRegion(const KernelRegion& cpu_kernel, const CudaPointer<CudaShell>& shell,
                                       const CudaPointer<CudaTriangulatedShell>& triangulation)
        : mShell(shell),
          mTriangulation(triangulation)
    {
        Eigen::MatrixXd A;
        zju_mat_to_eigen(cpu_kernel.A, A);
        this->A.assign(A);

        Eigen::MatrixXi lines;
        zju_mat_to_eigen(cpu_kernel.mLines, lines);
        mLines.assign(lines);

        this->mPointType = cpu_kernel.mPointType;
        this->mInvalidRegionType = cpu_kernel.mInvalidRegionType;
        this->mClockwise = cpu_kernel.mClockwise;

        for(size_t sample : cpu_kernel.mInnerSamples)
        {
            mInnerSamples.push_back(sample);
        }

        for(size_t sample : cpu_kernel.mOuterSamples)
        {
            mOuterSamples.push_back(sample);
        }
    }

    __device__ bool CudaKernelRegion::contains(const Eigen::Vector2d &point) const
    {
        Eigen::Vector3d homo;
        homo[0] = point[0];
        homo[1] = point[1];
        homo[2] = 1;

        const Eigen::MatrixXd& A = *this->A;
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

    __device__ bool CudaKernelRegion::isInvalidRegion(const Eigen::Vector2d& point) const
    {
        auto& lines = *mLines;
        auto& vertices = *mTriangulation->vertices;
        auto& innerSamples = mInnerSamples;
        auto& outerSamples = mOuterSamples;

        for(int i = 0; i < lines.cols(); i++)
        {
            for(int j = 0; j < innerSamples.size(); j++)
            {
                size_t sample = innerSamples[i];
                Eigen::Vector3d bary;

                if(barycentric_2D(vertices.col(lines(0, i)), vertices.col(lines(1, i)), point,
                                  mShell->innerShell->col(sample), bary))
                {
                    //the point is inside the tetrahedron
                    double f0 = mTriangulation->getFValue(lines(0, i));
                    double f1 = mTriangulation->getFValue(lines(1, i));
                    double f2 = mTriangulation->getFValue(mPointType);

                    double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2];
                    if(f > 0)
                    {
                        return true;
                    }
                }
            }

            for(int i = 0; i < outerSamples.size(); i++)
            {
                size_t sample = outerSamples[i];
                Eigen::Vector3d bary;
                if(barycentric_2D(vertices.col(lines(0, i)), vertices.col(lines(1, i)), point,
                                  mShell->outerShell->col(sample), bary))
                {
                    //the point is inside the tetrahedron
                    double f0 = mTriangulation->getFValue(lines(0, i));
                    double f1 = mTriangulation->getFValue(lines(1, i));
                    double f2 = mTriangulation->getFValue(mPointType);

                    double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2];
                    if(f < 0)
                    {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    __device__ bool CudaKernelRegion::barycentric(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                                                  const Eigen::Vector3d& p, Eigen::Vector3d& bary) const
    {
        const Eigen::Vector3d u = b - a;
        const Eigen::Vector3d v = c - a;
        const Eigen::Vector3d w = p - a;
        const Eigen::Vector3d vw = v.cross(w);
        const Eigen::Vector3d vu = v.cross(u);
        const Eigen::Vector3d uw = u.cross(w);
        const double denom = 1.0 / vu.norm();
        bary[1] = vu.dot(vw) >= 0 ? vw.norm() * denom : -vw.norm() * denom;
        bary[2] = uw.dot(vu) <= 0 ? uw.norm() * denom : -uw.norm() * denom;
        bary[0] = 1 - bary[1] - bary[2];

        bool result = (bary[0] > 0 || fabs(bary[0]) < 1e-6)
                && (bary[1] >0 || fabs(bary[1]) < 1e-6)
                && (bary[2] > 0 || fabs(bary[2]) < 1e-6);
        return result;
    }

    __device__ bool CudaKernelRegion::barycentric_2D(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c,
                                                     const Eigen::Vector2d& p, Eigen::Vector3d& bary) const
    {
        Eigen::Vector3d a_3d, b_3d, c_3d, p_3d;
        a_3d[0] = a[0];
        a_3d[1] = a[1];
        b_3d[0] = b[0];
        b_3d[1] = b[1];
        c_3d[0] = c[0];
        c_3d[1] = c[1];
        p_3d[0] = p[0];
        p_3d[1] = p[1];

        return barycentric(a_3d, b_3d, c_3d, p_3d, bary);
    }
}
