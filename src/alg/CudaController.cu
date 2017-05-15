#include "Common.h"
#include "CudaController.h"
#include "CudaControllerImpl.h"
#include "CudaSamplingTree.h"
#include "CudaTriangulatedShell.h"
#include "TriangulatedShell.h"
#include "KernelRegion.h"
#include <wkylib/Geometry/Util.h>
#include <limits>

using namespace WKYLIB::Geometry;
using namespace WKYLIB::Cuda;

namespace SBV
{
    static const int THREADS_PER_BLOCK = 512;

    __device__ double computeError(const Eigen::Matrix3d& Q, const Eigen::Vector2d& point)
    {
        Eigen::Vector3d p;
        p[0] = point[0];
        p[1] = point[1];
        p[2] = 1;
        Eigen::Vector3d temp = Q * p;
        return p[0] * temp[0] + p[1] * temp[1] + p[2] * temp[2];
    }

    __device__ void reduct(double cache_error[], Eigen::Vector2d cache_pos[])
    {
        int i = blockDim.x / 2;
        while(i)
        {
            if(threadIdx.x < i)
            {
                if(cache_error[threadIdx.x + i] < cache_error[threadIdx.x])
                {
                     cache_error[threadIdx.x] = cache_error[threadIdx.x + i];
                     cache_pos[threadIdx.x] = cache_pos[threadIdx.x + i];
                }
            }
            __syncthreads();
            i /= 2;
        }
    }

    //kernel functions
    __global__ void kernel_find_collapse_pos_boundary(
            CudaKernelRegion* kernel,
            bool* isInner,
            Eigen::Matrix3d* Q1,
            Eigen::Matrix3d* Q2,
            Eigen::Vector2d* out_pos,
            double* out_error)
    {
        __shared__ double cache_error[THREADS_PER_BLOCK];
        __shared__ Eigen::Vector2d cache_pos[THREADS_PER_BLOCK];

        const CudaVector<size_t>* samples_id_ptr = nullptr;
        const Eigen::Map<Eigen::MatrixXd>* samples_ptr = nullptr;
        if(*isInner)
        {
            samples_id_ptr = &kernel->getInnerSamples();
            samples_ptr = kernel->getShell()->innerShell.get();
        }
        else
        {
            samples_id_ptr = &kernel->getOuterSamples();
            samples_ptr = kernel->getShell()->outerShell.get();
        }
        const CudaVector<size_t>& samples = *samples_id_ptr;

        int tid = threadIdx.x;
        cache_error[threadIdx.x] = std::numeric_limits<double>::max();
        while(tid < samples.size())
        {
            int sampleId = samples[tid];
            const Eigen::Vector2d& candidate = samples_ptr->col(sampleId);
            if(kernel->contains(candidate))
            {
                double error = computeError(*Q1, candidate) + computeError(*Q2, candidate);;
                if(error < cache_error[threadIdx.x])
                cache_error[threadIdx.x] = error;
                cache_pos[threadIdx.x] = candidate;
            }

            tid += blockDim.x;
        }      
        __syncthreads();

        reduct(cache_error, cache_pos);

        if(threadIdx.x == 0)
        {
            if(cache_error[threadIdx.x] < std::numeric_limits<double>::max())
            {
                *out_error = cache_error[threadIdx.x];
                *out_pos = cache_pos[threadIdx.x];
            }
        }
    }

    __global__ void kernel_find_collapse_pos_zero_set(
            CudaKernelRegion* kernel,
            double* xmin,
            double* ymin,
            int* xCount,
            int* yCount,
            double* sampleRadius,
            Eigen::Matrix3d* Q1,
            Eigen::Matrix3d* Q2,
            Eigen::Vector2d* out_pos,
            double* out_error)
    {
        __shared__ double cache_error[THREADS_PER_BLOCK];
        __shared__ Eigen::Vector2d cache_pos[THREADS_PER_BLOCK];

        cache_error[threadIdx.x] = std::numeric_limits<double>::max();
        int tid = threadIdx.x;
        int total = *xCount * (*yCount);
        while(tid < total)
        {
            Eigen::Vector2d point;
            int y = tid / (*xCount);
            int x = tid - y * (*xCount);
            point[0] = *xmin + *sampleRadius * x;
            point[1] = *ymin + *sampleRadius * y;

            if(kernel->contains(point))
            {
                double error = computeError(*Q1, point) + computeError(*Q2, point);
                if(error < cache_error[threadIdx.x])
                cache_error[threadIdx.x] = error;
                cache_pos[threadIdx.x] = point;
            }
            tid += blockDim.x;
        }
        __syncthreads();

        reduct(cache_error, cache_pos);

        if(threadIdx.x == 0)
        {
            if(cache_error[threadIdx.x] < std::numeric_limits<double>::max())
            {
                *out_error = cache_error[threadIdx.x];
                *out_pos = cache_pos[threadIdx.x];
            }
        }
    }

    ///////////////////////////////////////////////////////////////////
    CudaController::~CudaController()
    {
        if(impl)
        {
            delete impl;
        }
    }

    void CudaController::build(const matrixr_t &innerShell, const matrixr_t &outerShell, const TriangulatedShell &triangulation)
    {
        if(impl)
        {
            delete impl;
        }
        this->impl = new CudaControllerImpl();
        this->impl->build(innerShell, outerShell, triangulation);
    }

    void CudaController::sample(double xmin, double xmax, double ymin, double ymax, double sampleRadius,
                                std::vector<matrixr_t> &output_samples)
    {
        this->impl->sample(xmin, xmax, ymin, ymax, sampleRadius, output_samples);
    }

    void CudaController::buildKernelRegion(const KernelRegion &kernel)
    {
        this->impl->buildKernelRegion(kernel);
    }

    bool CudaController::findCollapsePos_Boundary(bool isInner, const Eigen::Matrix3d& Q1, const Eigen::Matrix3d& Q2,\
                                                  matrixr_t& position, double& out_error)
    {
        return this->impl->findCollapsePos_Boundary(isInner, Q1, Q2, position, out_error);
    }

    bool CudaController::findCollapsePos_ZeroSet(const Eigen::Matrix3d& Q1, const Eigen::Matrix3d& Q2, double sampleRadius,
                                                 matrixr_t& position, double& out_error)
    {
        return this->impl->findCollapsePos_ZeroSet(Q1, Q2, sampleRadius, position, out_error);
    }

    //////////////////////////////////////////////////////////////////////////////////
    CudaControllerImpl::CudaControllerImpl()
    {
        cudaGetDeviceCount(&mDeviceCount);
        for(int i = 0; i < mDeviceCount; i++)
        {
            cudaGetDeviceProperties(&mDeviceProperty, i);
        }
    }

    void CudaControllerImpl::build(const matrixr_t &innerShell, const matrixr_t &outerShell, const TriangulatedShell &triangulation)
    {
        mShell.assign(CudaShell());

        Eigen::MatrixXd eigen_innerShell;
        Eigen::MatrixXd eigen_outerShell;
        zju_mat_to_eigen(innerShell, eigen_innerShell);
        zju_mat_to_eigen(outerShell, eigen_outerShell);
        mShell->innerShell.assign(eigen_innerShell);
        mShell->outerShell.assign(eigen_outerShell);

        castTriangulation(triangulation, mTriangulation);
    }

    void CudaControllerImpl::castTriangulation(const TriangulatedShell &triangulation, CudaPointer<CudaTriangulatedShell> &cuda_triangulation)
    {
        cuda_triangulation.assign(CudaTriangulatedShell());

        cuda_triangulation->vertices.reserve(triangulation.vertices.size(2));
        for(int i = 0; i < triangulation.vertices.size(2); i++)
        {
            Eigen::Vector2d vert;
            vert[0] = triangulation.vertices(0, i);
            vert[1] = triangulation.vertices(1, i);
            cuda_triangulation->vertices.push_back(vert);
        }

        cuda_triangulation->triangles.reserve(triangulation.triangles.size(2));
        for(int i = 0; i < triangulation.triangles.size(2); i++)
        {
            Eigen::Vector3i tri;
            tri[0] = triangulation.triangles(0, i);
            tri[1] = triangulation.triangles(1, i);
            tri[2] = triangulation.triangles(2, i);
            cuda_triangulation->triangles.push_back(tri);
        }

        cuda_triangulation->vertType.assign(triangulation.vertType);
    }

    void CudaControllerImpl::buildKernelRegion(const KernelRegion &kernel)
    {
        mKernel.assign(CudaKernelRegion(kernel, mShell, mTriangulation));
    }

    void CudaControllerImpl::sample(double xmin, double xmax, double ymin, double ymax, double sampleRadius,
                                    std::vector<matrixr_t> &output_samples)
    {
        CudaSamplingTree tree(mKernel, xmax, xmin, ymax, ymin, sampleRadius);

        output_samples.reserve(tree.getSamples().size());
        for(int i = 0; i < tree.getSamples().size(); i++)
        {
            matrixr_t mat;
            eigen_to_zju_vector2d(tree.getSamples()[i], mat);
            output_samples.push_back(mat);
        }
    }

    bool CudaControllerImpl::findCollapsePos_Boundary(bool isInner, const Eigen::Matrix3d& Q1, const Eigen::Matrix3d& Q2,
                                  matrixr_t& position, double& out_error)
    {
        //transfer data to gpu
        CudaPointer<bool> gpu_isInner(isInner);
        CudaPointer<double> gpu_outError(out_error);
        CudaPointer<Eigen::Matrix3d> gpu_Q1(Q1);
        CudaPointer<Eigen::Matrix3d> gpu_Q2(Q2);
        CudaPointer<Eigen::Vector2d> gpu_pos;
        gpu_pos.assign(Eigen::Vector2d());

        kernel_find_collapse_pos_boundary <<<1, THREADS_PER_BLOCK>>> (mKernel.get(), gpu_isInner.get(),
                                                                      gpu_Q1.get(), gpu_Q2.get(),
                                                                      gpu_pos.get(), gpu_outError.get());
        cudaDeviceSynchronize();

        if(*gpu_outError == std::numeric_limits<double>::max())
            return false;    //no position is inside valid region

        out_error = *gpu_outError;
        eigen_to_zju_vector2d(*gpu_pos, position);
        return true;
    }

    bool CudaControllerImpl::findCollapsePos_ZeroSet(const Eigen::Matrix3d& Q1, const Eigen::Matrix3d& Q2,
                                                     double sampleRadius,
                                                     matrixr_t& position,
                                                     double& out_error)
    {
        CudaPointer<double> gpu_outError(out_error);
        CudaPointer<Eigen::Matrix3d> gpu_Q1(Q1);
        CudaPointer<Eigen::Matrix3d> gpu_Q2(Q2);
        CudaPointer<Eigen::Vector2d> gpu_pos;
        gpu_pos.assign(Eigen::Vector2d());

        //find AABB of the one ring area
        double xmin = std::numeric_limits<double>::max();
        double xmax = -std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();
        double ymax = -std::numeric_limits<double>::max();

        for(int i = 0; i < mKernel->mLines.size(); i++)
        {
            const Eigen::Vector2d& a = mTriangulation->vertices[mKernel->mLines[i][0]];
            const Eigen::Vector2d& b = mTriangulation->vertices[mKernel->mLines[i][1]];
            xmax = a[0] > xmax ? a[0] : xmax;
            xmin = a[0] < xmin ? a[0] : xmin;
            ymax = a[1] > ymax ? a[1] : ymax;
            ymin = a[1] < ymin ? a[1] : ymin;
            xmax = b[0] > xmax ? b[0] : xmax;
            xmin = b[0] < xmin ? b[0] : xmin;
            ymax = b[1] > ymax ? b[1] : ymax;
            ymin = b[1] < ymin ? b[1] : ymin;
        }

        int xCount = (xmax - xmin) / sampleRadius;
        int yCount = (ymax - ymin) / sampleRadius;

        CudaPointer<int> gpu_xCount(xCount);
        CudaPointer<int> gpu_yCount(yCount);
        CudaPointer<double> gpu_xmin(xmin);
        CudaPointer<double> gpu_ymin(ymin);
        CudaPointer<double> gpu_sampleRadius(sampleRadius);

        kernel_find_collapse_pos_zero_set <<<1, THREADS_PER_BLOCK>>> (mKernel.get(), gpu_xmin.get(), gpu_ymin.get(), gpu_xCount.get(), gpu_yCount.get(),
                                          gpu_sampleRadius.get(), gpu_Q1.get(), gpu_Q2.get(), gpu_pos.get(), gpu_outError.get());
        cudaDeviceSynchronize();

        if(*gpu_outError == std::numeric_limits<double>::max())
            return false;    //no position is inside valid region

        out_error = *gpu_outError;
        eigen_to_zju_vector2d(*gpu_pos, position);
        return true;
    }
}
