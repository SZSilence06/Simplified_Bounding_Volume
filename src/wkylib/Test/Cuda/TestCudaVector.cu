#include <wkylib/Test/catch.hpp>
#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>

using namespace WKYLIB::Cuda;

TEST_CASE("CudaVector copy", "[CudaVector]") {
    CudaVector<int> objects;
    for(int i = 0; i < 10; i++)
    {
        objects.push_back(i);
    }

    REQUIRE(objects.size() == 10);
    for(int i = 0; i < objects.size(); i++)
    {
        CHECK(objects[i] == i);
    }

    CudaVector<int> anotherVector = objects;
    REQUIRE(anotherVector.size() == 10);
    for(int i = 0; i < anotherVector.size(); i++)
    {
        CHECK(anotherVector[i] == i);
    }
}

__global__ void kernel_test_vector_manipulation(CudaPointer<CudaVector<int>> vector)
{
    CudaVector<int>& v = *vector;
    for(int i = 0; i < v.size(); i++)
    {
        v[i] *= 2;
    }
}

TEST_CASE("CudaVector manipulation with gpu kernel", "[CudaVector]") {
    CudaVector<int> objects;
    for(int i = 0; i < 10; i++)
    {
        objects.push_back(i);
    }

    REQUIRE(objects.size() == 10);
    for(int i = 0; i < objects.size(); i++)
    {
        CHECK(objects[i] == i);
    }

    CudaPointer<CudaVector<int>> gpu_vector;
    gpu_vector.assign(objects);

    REQUIRE(gpu_vector->size() == 10);
    for(int i = 0; i < gpu_vector->size(); i++)
    {
        CudaVector<int>& v = *gpu_vector;
        CHECK(v[i] == i);
    }

    kernel_test_vector_manipulation <<< 1, 1 >>> (gpu_vector);
    cudaDeviceSynchronize();

    CudaVector<int>& v = *gpu_vector;
    for(int i = 0; i <v.size(); i++)
    {
        CHECK(v[i] == 2 * i);
    }
}
