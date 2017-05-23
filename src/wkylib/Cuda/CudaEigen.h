#ifndef WKY_CUDA_EIGEN_H
#define WKY_CUDA_EIGEN_H

#include <wkylib/External/eigen3.3/Eigen/Dense>
#include "CudaAllocator.h"
#include "CudaPointer.h"

namespace WKYLIB
{
    namespace Cuda
    {
        template<>
        class CudaAllocator<Eigen::MatrixXd>
        {
        public:
            CudaAllocator() = default;
            ~CudaAllocator() = default;

            static Eigen::MatrixXd* allocate()
            {
                Eigen::MatrixXd* result = nullptr;
                cudaMallocManaged(&result, sizeof(Eigen::Map<Eigen::MatrixXd>));
                return result;
            }

            static Eigen::MatrixXd* allocate(int count)
            {
                Eigen::MatrixXd* result = nullptr;
                cudaMallocManaged(&result, sizeof(Eigen::Map<Eigen::MatrixXd>) * count);
                return result;
            }

            static void construct(Eigen::MatrixXd* object, const Eigen::MatrixXd& other)
            {
                const double* cpu_data = other.data();
                double* gpu_data = nullptr;
                cudaMallocManaged(&gpu_data, sizeof(double) * other.cols() * other.rows());
                memcpy(gpu_data, cpu_data, sizeof(double) * other.cols() * other.rows());
                new (object) Eigen::Map<Eigen::MatrixXd>(gpu_data, other.rows(), other.cols());
            }

            template<class... Args >
            static void construct(Eigen::MatrixXd* object, Args&& ... args)
            {
                new (object) Eigen::MatrixXd(std::forward<Args>(args)...);
            }

            static void destroy(Eigen::MatrixXd* object)
            {
                cudaFree(object->data());
            }

            static void deallocate(Eigen::MatrixXd* memory)
            {
                cudaFree(memory);
            }
        };

        template<>
        class CudaAllocator<Eigen::MatrixXi>
        {
        public:
            CudaAllocator() = default;
            ~CudaAllocator() = default;

            static Eigen::MatrixXi* allocate()
            {
                Eigen::MatrixXi* result = nullptr;
                cudaMallocManaged(&result, sizeof(Eigen::Map<Eigen::MatrixXi>));
                return result;
            }

            static Eigen::MatrixXi* allocate(int count)
            {
                Eigen::MatrixXi* result = nullptr;
                cudaMallocManaged(&result, sizeof(Eigen::Map<Eigen::MatrixXi>) * count);
                return result;
            }

            static void construct(Eigen::MatrixXi* object, const Eigen::MatrixXi& other)
            {
                const int* cpu_data = other.data();
                int* gpu_data = nullptr;
                cudaMallocManaged(&gpu_data, sizeof(int) * other.cols() * other.rows());
                memcpy(gpu_data, cpu_data, sizeof(int) * other.cols() * other.rows());
                new (object) Eigen::Map<Eigen::MatrixXi>(gpu_data, other.rows(), other.cols());
            }

            template<class... Args >
            static void construct(Eigen::MatrixXi* object, Args&& ... args)
            {
                new (object) Eigen::MatrixXi(std::forward<Args>(args)...);
            }

            static void destroy(Eigen::MatrixXi* object)
            {
                cudaFree(object->data());
            }

            static void deallocate(Eigen::MatrixXi* memory)
            {
                cudaFree(memory);
            }
        };

        template<>
        class CudaPointer<Eigen::MatrixXd>
        {
        public:
            CudaPointer()
            {

            }

            CudaPointer(const Eigen::MatrixXd& obj) : CudaPointer()
            {
                assign(obj);
            }

            CudaPointer(const CudaPointer<Eigen::MatrixXd>& another)
            {
                this->refCount = another.refCount;
                this->pointer = another.pointer;
                if(this->refCount)
                {
                    this->refCount->addRef();
                }
            }

            CudaPointer(CudaPointer<Eigen::MatrixXd>&& rhs)
            {
                this->refCount = rhs.refCount;
                this->pointer = rhs.pointer;
                if(this->refCount)
                {
                    this->refCount->addRef();
                }
            }

            CudaPointer& operator= (const CudaPointer<Eigen::MatrixXd>& another)
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

            void assign(const Eigen::MatrixXd& obj)
            {
                if(this->refCount)
                {
                    this->refCount->decRef();
                }
                this->refCount = new RefCount<Eigen::MatrixXd>(obj);
                this->pointer = (Eigen::Map<Eigen::MatrixXd>*)refCount->get();
            }

            Eigen::Map<Eigen::MatrixXd>& __host__ __device__ operator*()
            {
                return *this->pointer;
            }

            const Eigen::Map<Eigen::MatrixXd>& __host__ __device__ operator*() const
            {
                return *this->pointer;
            }

            __host__ __device__ Eigen::Map<Eigen::MatrixXd>* operator->()
            {
                return this->pointer;
            }

            __host__ __device__ const Eigen::Map<Eigen::MatrixXd>* operator->() const
            {
                return this->pointer;
            }

            __host__ __device__ Eigen::Map<Eigen::MatrixXd>* get()
            {
                return this->pointer;
            }

            __host__ __device__ const Eigen::Map<Eigen::MatrixXd>* get() const
            {
                return this->pointer;
            }

        private:
            Eigen::Map<Eigen::MatrixXd>* pointer = nullptr;
            RefCount<Eigen::MatrixXd>* refCount = nullptr;
        };

        template<>
        class CudaPointer<Eigen::MatrixXi>
        {
        public:
            CudaPointer()
            {

            }

            CudaPointer(const Eigen::MatrixXi& obj) : CudaPointer()
            {
                assign(obj);
            }

            CudaPointer(const CudaPointer<Eigen::MatrixXi>& another)
            {
                this->refCount = another.refCount;
                this->pointer = another.pointer;
                if(this->refCount)
                {
                    this->refCount->addRef();
                }
            }

            CudaPointer(CudaPointer<Eigen::MatrixXi>&& rhs)
            {
                this->refCount = rhs.refCount;
                this->pointer = rhs.pointer;
                if(this->refCount)
                {
                    this->refCount->addRef();
                }
            }

            CudaPointer& operator= (const CudaPointer<Eigen::MatrixXi>& another)
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

            void assign(const Eigen::MatrixXi& obj)
            {
                if(this->refCount)
                {
                    this->refCount->decRef();
                }
                this->refCount = new RefCount<Eigen::MatrixXi>(obj);
                this->pointer = (Eigen::Map<Eigen::MatrixXi>*)refCount->get();
            }

            Eigen::Map<Eigen::MatrixXi>& __host__ __device__ operator*()
            {
                return *this->pointer;
            }

            const Eigen::Map<Eigen::MatrixXi>& __host__ __device__ operator*() const
            {
                return *this->pointer;
            }

            __host__ __device__ Eigen::Map<Eigen::MatrixXi>* operator->()
            {
                return this->pointer;
            }

            __host__ __device__ const Eigen::Map<Eigen::MatrixXi>* operator->() const
            {
                return this->pointer;
            }

            __host__ __device__ Eigen::Map<Eigen::MatrixXi>* get()
            {
                return this->pointer;
            }

            __host__ __device__ const Eigen::Map<Eigen::MatrixXi>* get() const
            {
                return this->pointer;
            }

        private:
            Eigen::Map<Eigen::MatrixXi>* pointer = nullptr;
            RefCount<Eigen::MatrixXi>* refCount = nullptr;
        };
    }
}

#endif
