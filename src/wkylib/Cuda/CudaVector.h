#ifndef WKY_CUDA_VECTOR_H
#define WKY_CUDA_VECTOR_H

#include "CudaAllocator.h"
#include <vector>
#include <cassert>

namespace WKYLIB
{
    namespace Cuda
    {
        template<class T>
        class CudaVector
        {
        public:
            using iterator = T*;
            using const_iterator = const T*;
            using pointer = T*;
            using reference = T*;
            using allocator = CudaAllocator<T>;

        public:
            __host__ __device__ CudaVector()
            {
            }

            CudaVector(const CudaVector<T>& other)
            {
                mSize = other.mSize;
                mCapacity = other.mCapacity;
                mElements = allocator::allocate(mCapacity);
                for(int i = 0; i < mSize; i++)
                {
                    allocator::construct(mElements + i, other[i]);
                }
            }

            CudaVector(const std::vector<T>& other)
            {
                assign(other);
            }

            CudaVector(CudaVector<T>&& rhs)
            {
                std::swap(mSize, rhs.mSize);
                std::swap(mCapacity, rhs.mCapacity);
                std::swap(mElements, rhs.mElements);
            }

            ~CudaVector()
            {
                destroyAll();
            }

            CudaVector& operator =(const CudaVector<T>& other)
            {
                destroyAll();

                mSize = other.mSize;
                mCapacity = other.mCapacity;
                mElements = allocator::allocate(mCapacity);
                for(int i = 0; i < mSize; i++)
                {
                    allocator::construct(mElements + i, other[i]);
                }
                return *this;
            }

            CudaVector& operator =(CudaVector<T>&& rhs)
            {
                std::swap(mSize, rhs.mSize);
                std::swap(mCapacity, rhs.mCapacity);
                std::swap(mElements, rhs.mElements);
                return *this;
            }

            iterator begin()
            {
                return mElements;
            }

            iterator end()
            {
                return mElements + mSize + 1;
            }

            const_iterator cbegin() const
            {
                return mElements;
            }

            const_iterator cend() const
            {
                return mElements + mSize + 1;
            }

            void push_back(const T& object)
            {
                if(mSize >= mCapacity)
                {
                    expandCapacity();
                }
                allocator::construct(mElements + mSize, object);
                mSize++;
            }

            __device__ void push_back_on_device(const T& object)
            {
                assert(mSize < mCapacity);
                allocator::construct(mElements + mSize, object);
                mSize++;
            }

            void reserve(size_t size)
            {
                destroyAll();
                mCapacity = size;
                mElements = allocator::allocate(mCapacity);
            }

            void assign(const std::vector<T>& vector)
            {
                destroyAll();
                mSize = vector.size();
                mCapacity = vector.capacity();
                mElements = allocator::allocate(mCapacity);
                for(int i = 0; i < vector.size(); i++)
                {
                    allocator::construct(mElements + i, vector[i]);
                }
            }

            __host__ __device__ size_t size() const { return mSize; }
            __host__ __device__ size_t capacity() const { return mCapacity; }

            __host__ __device__ T& operator[] (int index)
            {
                assert(index >= 0 && index < mSize);
                return mElements[index];
            }

            __host__ __device__ const T& operator[] (int index) const
            {
                assert(index >= 0 && index < mSize);
                return mElements[index];
            }

        private:
            void expandCapacity()
            {
                T* old = mElements;
                if(mCapacity == 0)
                {
                    mCapacity = 1;
                }
                else
                {
                    mCapacity *= 2;
                }
                cudaMallocManaged(&mElements, sizeof(T) * mCapacity);
                for(int i = 0; i < mSize; i++)
                {
                    allocator::construct(mElements + i, std::move(old[i]));
                    allocator::destroy(old + i);
                }
                allocator::deallocate(old);
            }

            void destroyAll()
            {
                if(mElements)
                {
                    for(int i = 0; i < mSize; i++)
                    {
                        allocator::destroy(mElements + i);
                    }
                    allocator::deallocate(mElements);
                    mElements = nullptr;
                    mSize = mCapacity = 0;
                }
            }

        private:
            T* mElements = nullptr;
            size_t mSize = 0;
            size_t mCapacity = 0;
        };

        //template specialization
        template<class T>
        class CudaAllocator<CudaVector<T>>
        {
        public:
            CudaAllocator() = default;
            ~CudaAllocator() = default;

            static CudaVector<T>* allocate()
            {
                CudaVector<T>* result = nullptr;
                cudaMallocManaged(&result, sizeof(CudaVector<T>));
                return result;
            }

            static CudaVector<T>* allocate(int count)
            {
                CudaVector<T>* result = nullptr;
                if(count)
                    cudaMallocManaged(&result, sizeof(CudaVector<T>) * count);
                return result;
            }

            template<class... Args >
            static void construct(CudaVector<T>* object, Args&& ... args)
            {
                new (object) CudaVector<T>(std::forward<Args>(args)...);
            }

            static void destroy(CudaVector<T>* object)
            {
                object->~CudaVector();
            }

            static void deallocate(CudaVector<T>* memory)
            {
                if(memory)
                    cudaFree(memory);
            }
        };
    }
}

#endif

