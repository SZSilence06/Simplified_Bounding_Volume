#ifndef WKY_CUDA_POINTER_H
#define WKY_CUDA_POINTER_H

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
            RefCount(const T& object)
            {
                cudaMallocManaged(&this->object, sizeof(T));
                new(this->object) T(object);
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
                this->object->~T();
                cudaFree(this->object);
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
                    destroy();
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
            __host__ __device__ CudaPointer()
            {

            }

            __host__ __device__ CudaPointer(const T& obj) : CudaPointer()
            {
                assign(obj);
            }

            __host__ __device__ CudaPointer(const CudaPointer<T>& another)
            {
                this->refCount = another.refCount;
                this->pointer = another.pointer;
                if(this->refCount)
                {
                    this->refCount->addRef();
                }
            }

            __host__ __device__ CudaPointer(CudaPointer<T>&& rhs)
            {
                this->refCount = rhs.refCount;
                this->pointer = rhs.pointer;
                if(this->refCount)
                {
                    this->refCount->addRef();
                }
            }

            __host__ __device__ CudaPointer& operator= (const CudaPointer<T>& another)
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

            __host__ __device__ ~CudaPointer()
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
            void destroy()
            {
                if(this->pointer)
                {
                    this->pointer->~T();
                    cudaFree(this->pointer);
                    delete refCount;
                }
            }

        private:
            T* pointer = nullptr;
            RefCount<T>* refCount = nullptr;
        };
    }
}

#endif
