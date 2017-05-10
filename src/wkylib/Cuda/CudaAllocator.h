#ifndef WKY_CUDA_ALLOCATOR_H
#define WKY_CUDA_ALLOCATOR_H

namespace WKYLIB
{
    namespace Cuda
    {
        template<class T>
        class CudaAllocator
        {
        public:
            CudaAllocator() = default;
            ~CudaAllocator() = default;

            static T* allocate()
            {
                T* result = nullptr;
                cudaMallocManaged(&result, sizeof(T));
                return result;
            }

            static T* allocate(int count)
            {
                T* result = nullptr;
                if(count)
                    cudaMallocManaged(&result, sizeof(T) * count);
                return result;
            }

            template<class... Args >
            __host__ __device__ static void construct(T* object, Args&& ... args)
            {
                new (object) T(std::forward<Args>(args)...);
            }

            static void destroy(T* object)
            {
                object->~T();
            }

            static void deallocate(T* memory)
            {
                if(memory)
                    cudaFree(memory);
            }
        };
    }
}

#endif
