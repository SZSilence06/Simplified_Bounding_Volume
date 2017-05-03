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

class TestCudaVectorObject
{
public:
    CudaVector<int> data;
};

__global__ void kernel_test_vector_object_manipulation(CudaPointer<TestCudaVectorObject> object)
{
    CudaVector<int>& v = object->data;
    for(int i = 0; i < v.size(); i++)
    {
        v[i] *= 2;
    }
}

TEST_CASE("CudaVector-contained object manipulation with gpu kernel", "[CudaVector]") {
    TestCudaVectorObject object;
    for(int i = 0; i < 10; i++)
    {
        object.data.push_back(i);
    }

    REQUIRE(object.data.size() == 10);
    for(int i = 0; i < object.data.size(); i++)
    {
        CHECK(object.data[i] == i);
    }

    CudaPointer<TestCudaVectorObject> gpu_object;
    gpu_object.assign(object);

    REQUIRE(gpu_object->data.size() == 10);
    for(int i = 0; i < gpu_object->data.size(); i++)
    {
        CudaVector<int>& v = gpu_object->data;
        CHECK(v[i] == i);
    }

    kernel_test_vector_object_manipulation <<< 1, 1 >>> (gpu_object);
    cudaDeviceSynchronize();

    CudaVector<int>& v = gpu_object->data;
    for(int i = 0; i < v.size(); i++)
    {
        CHECK(v[i] == 2 * i);
    }
}
