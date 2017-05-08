#ifndef WKY_CUDA_EIGEN_H
#define WKY_CUDA_EIGEN_H

#include "Common.h"
#include "eigen3.3/Eigen/Dense"
#include <wkylib/Cuda/CudaAllocator.h>

template<>
class WKYLIB::Cuda::CudaAllocator<Eigen::MatrixXd>
{
public:
    CudaAllocator() = default;
    ~CudaAllocator() = default;

    static Eigen::MatrixXd* allocate()
    {
        Eigen::MatrixXd* result = nullptr;
        cudaMallocManaged(&result, sizeof(Eigen::MatrixXd));
        return result;
    }

    static Eigen::MatrixXd* allocate(int count)
    {
        Eigen::MatrixXd* result = nullptr;
        cudaMallocManaged(&result, sizeof(Eigen::MatrixXd) * count);
        return result;
    }

    static void construct(Eigen::MatrixXd* object, const Eigen::MatrixXd& other)
    {
        const double* cpu_data = other.data();
        double* gpu_data = nullptr;
        cudaMallocManaged(&gpu_data, sizeof(double) * other.size());
        memcpy(gpu_data, cpu_data, sizeof(double) * other.size());
        new (object) Eigen::Map<Eigen::MatrixXd>(gpu_data, other.rows(), other.cols());
    }

    template<class... Args >
    static void construct(Eigen::MatrixXd* object, Args&& ... args)
    {
        new (object) Eigen::MatrixXd(std::forward<Args>(args)...);
    }

    static void destroy(Eigen::MatrixXd* object)
    {
        cudaFree(object->data());
    }

    static void deallocate(Eigen::MatrixXd* memory)
    {
        cudaFree(memory);
    }
};

template<>
class WKYLIB::Cuda::CudaAllocator<Eigen::MatrixXi>
{
public:
    CudaAllocator() = default;
    ~CudaAllocator() = default;

    static Eigen::MatrixXi* allocate()
    {
        Eigen::MatrixXi* result = nullptr;
        cudaMallocManaged(&result, sizeof(Eigen::MatrixXi));
        return result;
    }

    static Eigen::MatrixXi* allocate(int count)
    {
        Eigen::MatrixXi* result = nullptr;
        cudaMallocManaged(&result, sizeof(Eigen::MatrixXi) * count);
        return result;
    }

    static void construct(Eigen::MatrixXi* object, const Eigen::MatrixXi& other)
    {
        const int* cpu_data = other.data();
        int* gpu_data = nullptr;
        cudaMallocManaged(&gpu_data, sizeof(int) * other.size());
        memcpy(gpu_data, cpu_data, sizeof(int) * other.size());
        new (object) Eigen::Map<Eigen::MatrixXi>(gpu_data, other.rows(), other.cols());
    }

    template<class... Args >
    static void construct(Eigen::MatrixXi* object, Args&& ... args)
    {
        new (object) Eigen::MatrixXi(std::forward<Args>(args)...);
    }

    static void destroy(Eigen::MatrixXi* object)
    {
        cudaFree(object->data());
    }

    static void deallocate(Eigen::MatrixXi* memory)
    {
        cudaFree(memory);
    }
};

#endif
