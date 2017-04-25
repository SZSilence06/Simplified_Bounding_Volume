#ifndef WKY_CUDA_VECTOR_H
#define WKY_CUDA_VECTOR_H

#include <thrust/device_vector.h>

namespace WKYLIB
{
    namespace Cuda
    {
        template<class T>
        class CudaVector
        {
        private:
            thrust::device_vector<T> mVector;
            T* mElements;
            size_t mSize;

        public:
            void assign(const std::vector<T>& vector)
            {
                mVector = vector;
                mElements = thrust::raw_pointer_cast(&mVector[0]);
                mSize = mVector.size();
            }

            __device__ size_t size() const { return mSize; }
            __device__ T* getElements() const { return mElements; }
        };
    }
}

#endif
