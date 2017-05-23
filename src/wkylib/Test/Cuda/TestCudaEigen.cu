#include <wkylib/Test/catch.hpp>
#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>
#include <wkylib/External/eigen3.3/Eigen/Dense>
#include <wkylib/Cuda/CudaEigen.h>

using namespace WKYLIB::Cuda;

__global__ void kernel__check__matrixXd(CudaPointer<Eigen::MatrixXd> gpu_mat, CudaPointer<CudaVector<double>> result)
{
    for(int i = 0; i < gpu_mat->rows(); i++)
    {
        for(int j = 0; j < gpu_mat->cols(); j++)
            result->push_back_on_device((*gpu_mat)(i, j));
    }
}

TEST_CASE("Eigen copy", "[CudaEigen]") {
    Eigen::MatrixXd cpu_mat;
    cpu_mat.resize(8, 8);
    for(int i = 0; i < cpu_mat.rows(); i++)
    {
        for(int j = 0; j < cpu_mat.cols(); j++)
            cpu_mat(i, j) = 8 * i + j;
    }

    CudaPointer<Eigen::MatrixXd> gpu_mat;
    gpu_mat.assign(cpu_mat);

    for(int i = 0; i < gpu_mat->rows(); i++)
    {
        for(int j = 0; j < gpu_mat->cols(); j++)
            CHECK((*gpu_mat)(i, j) == 8 * i + j);
    }

    CudaPointer<CudaVector<double>> result;
    result.assign(CudaVector<double>());
    result->reserve(100);
    kernel__check__matrixXd <<<1, 1>>> (gpu_mat, result);
    cudaDeviceSynchronize();
    for(int i = 0; i < cpu_mat.cols() * cpu_mat.rows(); i++)
    {
        CHECK((*result)[i] == i);
    }
}

