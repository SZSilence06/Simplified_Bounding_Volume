#include <wkylib/Test/catch.hpp>
#include <wkylib/Cuda/CudaPointer.h>

using namespace WKYLIB::Cuda;

TEST_CASE( "CudaPointer initialization", "[CudaPointer]" ) {
    int test = 20;

    SECTION("initialize with value")
    {
        CudaPointer<int> p(test);
        REQUIRE(*p == test);
    }

    SECTION("initialize with assign()")
    {
        CudaPointer<int> p;
        p.assign(test);
        REQUIRE(*p == test);
    }
}

class TestObject
{
public:
    CudaPointer<int> data;
};

__global__ void kernel_test_cuda_manipulation(CudaPointer<int> p)
{
    *p = 2;
}

__global__ void kernel_test_cuda_manipulation_with_object(CudaPointer<TestObject> p)
{
    *p->data = 2;
}

TEST_CASE( "CudaPointer manipulation with gpu kernel", "[CudaPointer]" ) {
    int test = 20;
    TestObject testObject;
    testObject.data.assign(test);

    CudaPointer<int> p(test);
    kernel_test_cuda_manipulation <<<1, 1>>> (p);
    cudaDeviceSynchronize();
    REQUIRE(*p == 2);

    CudaPointer<TestObject> p2(testObject);
    kernel_test_cuda_manipulation_with_object <<<1, 1>>> (p2);
    cudaDeviceSynchronize();
    REQUIRE(*p2->data == 2);
}






