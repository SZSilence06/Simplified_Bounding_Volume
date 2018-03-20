#ifndef WKY_SBV_GPU_FMM_H
#define WKY_SBV_GPU_FMM_H

#include "FMM.h"
#include <thrust/complex.h>
#include <eigen3/Eigen/Dense>
#include <boost/math/special_functions/legendre.hpp>
#include <wkylib/Cuda/CudaVector.h>
#include <thrust/complex.h>

namespace SBV
{
    class GPU_FMM;

    namespace __FMM_Internal{
        __global__ void kernel_downwardPass(GPU_FMM* fmm, int* level);
    }

    class GPU_MultipoleExp
    {
    public:
        __host__ __device__ thrust::complex<double>& moment(int n, int m) { return _moment[n][m + MAX_ORDER]; }
        __host__ __device__ const thrust::complex<double>& moment(int n, int m) const { return _moment[n][m + MAX_ORDER]; }
    private:
        thrust::complex<double> _moment[MAX_ORDER+1][MAX_ORDER * 2 + 1];
    };

    class GPU_LocalExp
    {
    public:
        __host__ __device__ thrust::complex<double>& moment(int n, int m) { return _moment[n][m + MAX_ORDER]; }
        __host__ __device__ const thrust::complex<double>& moment(int n, int m) const { return _moment[n][m + MAX_ORDER]; }
    private:
        thrust::complex<double> _moment[MAX_ORDER+1][MAX_ORDER * 2 + 1];
    };

    struct GPU_Face
    {
        Eigen::Matrix3d triangle;
        double derivative;
    };

    struct GPU_Cell
    {
        Eigen::Vector3d centroid;
        GPU_MultipoleExp multipoleExp;
        GPU_LocalExp localExp;
        WKYLIB::Cuda::CudaVector<GPU_Face> faces;
        GPU_Cell* parent = nullptr;
        WKYLIB::Cuda::CudaVector<GPU_Cell*> children;
        WKYLIB::Cuda::CudaVector<GPU_Cell*> interList;
        WKYLIB::Cuda::CudaVector<GPU_Cell*> neighbours;
        int xIndex, yIndex, zIndex;
        int level;
        bool hasFace;
    };

    /**
     * @brief The GPU_FMM class
     *
     *        This class is used to perform FMM field computation on GPU. And This class can only be used in cuda kernels.
     *
     *        To use this class, one should construct a @see FMM class first, and then construct a GPU_FMM class from it using @see buildFromCPU.
     */
    class GPU_FMM {
    public:
        ~GPU_FMM();

        void buildFromCPU(const FMM &fmm);

        __host__ __device__ double getPotential(const Eigen::Vector3d& x);

    private:
        void CpuCell2Gpu(CellPtr cpu_cell, GPU_Cell* gpu_cell);
        void CpuCellVec2Gpu(const std::vector<CellPtr>& cpu_vec, WKYLIB::Cuda::CudaVector<GPU_Cell*> &gpu_vec);

        __host__ __device__ double GPU_legendre_p(int n, int m, double x);
        __host__ __device__ inline static void GPU_toSphericalCoordinate(const Eigen::Vector3d& x, double& r, double& theta, double& phi);
        __host__ __device__ inline thrust::complex<double> GPU_R(const Eigen::Vector3d& x, int m, int n);
        __host__ __device__ inline thrust::complex<double> GPU_S(const Eigen::Vector3d& x, int m, int n);
        __host__ __device__ void GPU_M2L(const GPU_MultipoleExp &inputMoment, const Eigen::Vector3d &xc, const Eigen::Vector3d &x0, GPU_LocalExp &result);
        __host__ __device__ void GPU_L2L(const GPU_LocalExp &inputMoment, const Eigen::Vector3d &x0, const Eigen::Vector3d &x1, GPU_LocalExp &result);
        __host__ __device__ void getCellIndex(const Eigen::Vector3d& center, size_t& xIndex, size_t& yIndex, size_t& zIndex, size_t level);
        __host__ __device__ double localEvaluate(const GPU_LocalExp& mul, const Eigen::Vector3d& x, const Eigen::Vector3d& x0);
        __host__ __device__ double directEvaluate(const GPU_Face& face, const Eigen::Vector3d& x);

    private:
        int mMaxLevel = 0;

        WKYLIB::Cuda::CudaVector<double> mStep;
        double mXMax, mXMin, mYMax, mYMin, mZMax, mZMin;

        double factorial[100];
        double factorial_reciprocal[100];
        double double_factorial[100];

        GPU_Cell***** mCells;

        friend class FMM;

        friend __global__ void __FMM_Internal::kernel_downwardPass(GPU_FMM* fmm, int* level);
    };
}

#endif
