#ifndef WKY_CUDA_POINTER_H
#define WKY_CUDA_POINTER_H

namespace WKYLIB
{
    namespace Cuda
    {
        template<class T>
        class CudaPointer
        {
        public:
            CudaPointer()
            {
                cudaMallocManaged(&this->pointer, sizeof(T));
            }

            CudaPointer(const T& obj) : CudaPointer()
            {
                *this->pointer = obj;
            }

            CudaPointer(const CudaPointer<T>& another) = delete;
            CudaPointer(CudaPointer<T>&& rhs) = delete;

            ~CudaPointer()
            {
                cudaFree(this->pointer);
            }

            CudaPointer& operator= (const T& obj)
            {
                *this->pointer = obj;
            }

            CudaPointer& operator= (const CudaPointer<T>& another) = delete;

            T* operator->()
            {
                return this->pointer;
            }

            const T* operator->() const
            {
                return this->pointer;
            }

            T* get()
            {
                return this->pointer;
            }

            const T* get() const
            {
                return this->pointer;
            }

        private:
            T* pointer = nullptr;
        };
    }
}

#endif
