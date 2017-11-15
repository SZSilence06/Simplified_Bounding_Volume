#ifndef WKY_CUDA_POINTER_H
#define WKY_CUDA_POINTER_H

#include "CudaAllocator.h"

namespace WKYLIB
{
    namespace Cuda
    {
        template<class T>
        class CudaPointer;

        template<class T>
        class RefCount
        {
        private:
            using Allocator = CudaAllocator<T>;

            RefCount(const T& object)
            {
                this->object = Allocator::allocate();
                Allocator::construct(this->object, object);
                this->ref = 1;
            }

            RefCount(T&& object)
            {
                this->object = Allocator::allocate();
                Allocator::construct(this->object, std::forward<T&&>(object));
                this->ref = 1;
            }


            ~RefCount()
            {
                if(this->object)
                {
                    destroy();
                }
            }

            T* get() { return this->object;}

            void destroy()
            {
                Allocator::destroy(object);
                Allocator::deallocate(this->object);
                this->object = nullptr;
            }

            void addRef()
            {
                this->ref++;
            }

            void decRef()
            {
                this->ref--;
                if(this->ref == 0)
                {
                    delete this;
                }               
            }

        private:
            T* object = nullptr;
            int ref = 0;

            friend class CudaPointer<T>;
        };

        template<class T>
        class CudaPointer
        {
        public:
            CudaPointer()
            {

            }

            CudaPointer(const T& obj) : CudaPointer()
            {
                assign(obj);
            }

            CudaPointer(const CudaPointer<T>& another)
            {
                this->refCount = another.refCount;
                this->pointer = another.pointer;
                if(this->refCount)
                {
                    this->refCount->addRef();
                }
            }

            CudaPointer(CudaPointer<T>&& rhs)
            {
                this->refCount = rhs.refCount;
                this->pointer = rhs.pointer;
                if(this->refCount)
                {
                    this->refCount->addRef();
                }
            }

            CudaPointer& operator= (const CudaPointer<T>& another)
            {
                if(this == &another)
                {
                    return *this;
                }

                if(this->refCount)
                {
                    this->refCount->decRef();
                }

                this->refCount = another.refCount;
                this->pointer = another.pointer;
                if(this->refCount)
                {
                    this->refCount->addRef();
                }

                return *this;
            }

            ~CudaPointer()
            {
                if(this->refCount)
                {
                    this->refCount->decRef();
                }
            }

            void assign(const T& obj)
            {
                if(this->refCount)
                {
                    this->refCount->decRef();
                }
                this->refCount = new RefCount<T>(obj);
                this->pointer = refCount->get();
            }

            void assign(T&& obj)
            {
                if(this->refCount)
                {
                    this->refCount->decRef();
                }
                this->refCount = new RefCount<T>(std::forward<T&&>(obj));
                this->pointer = refCount->get();
            }

            T& __host__ __device__ operator*()
            {
                return *this->pointer;
            }

            const T& __host__ __device__ operator*() const
            {
                return *this->pointer;
            }

            __host__ __device__ T* operator->()
            {
                return this->pointer;
            }

            __host__ __device__ const T* operator->() const
            {
                return this->pointer;
            }

            __host__ __device__ T* get()
            {
                return this->pointer;
            }

            __host__ __device__ const T* get() const
            {
                return this->pointer;
            }

        private:
            T* pointer = nullptr;
            RefCount<T>* refCount = nullptr;
        };
    }
}

#endif
