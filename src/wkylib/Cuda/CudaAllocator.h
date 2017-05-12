#ifndef WKY_CUDA_ALLOCATOR_H
#define WKY_CUDA_ALLOCATOR_H

#ifndef __global__
#define __global__
#endif

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif

#ifndef __shared__
#define __shared__
#endif

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
