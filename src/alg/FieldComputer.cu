#include "FieldComputer.h"
#include <zjucad/matrix/matrix.h>
#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>
#include <eigen3/Eigen/Dense>
#include <iostream>

using namespace zjucad::matrix;
using namespace WKYLIB::Cuda;

namespace SBV {
    struct GPU_SamplePoint
    {
        Eigen::Vector3d position;
        Eigen::Vector3d normal;
        double value = 0;
        double derivative = 0;
        double size = 0;   //indicating the size of the triangle which the sample point lies in.
        Eigen::Matrix3d tri;      //indicating the triangle which the sample point lies in.
        Eigen::Matrix4d  transform;  //matrix for transforming to local
        Eigen::Matrix4d  invTransform;
    };

    const int THREADS_PER_BLOCK = 256;
    const int BLOCK_COUNT = 4;

    template<class T>
    __device__ static inline T my_min(const T& a, const T& b) {
        return a < b ? a : b;
    }

    __device__ static void integrateOverTriangle(const Eigen::Vector3d& x, const GPU_SamplePoint &point, double& I1, Eigen::Vector3d& Igrad)
    {
        Eigen::Matrix4d triangle;
        triangle.block<3, 3>(0, 0) = point.tri;
        triangle(3, 0) = triangle(3, 1) = triangle(3, 2) = 1;
        Eigen::Matrix4d localTriangle = point.transform * triangle;

        double l3 = localTriangle(0, 1);
        double u3 = localTriangle(0, 2);
        double v3 = localTriangle(1, 2);

        Eigen::Vector4d tempX;
        tempX.block<3, 1>(0, 0) = x;
        tempX[3] = 1;
        Eigen::Vector4d localX = point.transform * tempX;
        double u0 = localX[0];
        double v0 = localX[1];
        double w0 = localX[2];

        // edge lengths
        double l1 = sqrt((l3-u3) * (l3-u3) + v3*v3);
        double l2 = sqrt(u3*u3 + v3*v3);

        // threshold for small numbers
        double threshold = 1e-6 * my_min(my_min(l1,l2), l3);
        if(fabs(w0) < threshold)
            w0 = 0;

        // eq (3)
        Eigen::Vector3d sminus, splus;
        sminus[0] = -((l3-u3)*(l3-u0)+v3*v0)/l1;
        sminus[1] = -(u3*(u3-u0)+v3*(v3-v0))/l2;
        sminus[2] = -u0;
        splus[0] = ((u3-l3)*(u3-u0)+v3*(v3-v0))/l1;
        splus[1] = (u3*u0+v3*v0)/l2;
        splus[2] = l3-u0;

        // eq (4)
        Eigen::Vector3d t0;
        t0[0] = ((u3-l3)*v0+v3*(l3-u0))/l1;
        t0[1] = (v3*u0-u3*v0)/l2;
        t0[2] = v0;

        // eq (5)
        Eigen::Vector3d tplus, tminus;
        tplus[0] = sqrt((u3-u0)*(u3-u0) + (v3-v0)*(v3-v0));
        tplus[1] = sqrt(u0*u0 + v0*v0);
        tplus[2] = sqrt((l3-u0)*(l3-u0) + v0*v0);
        tminus[0] = tplus[2];
        tminus[1] = tplus[0];
        tminus[2] = tplus[1];

        // line 1, pp. 1450
        Eigen::Vector3d R0;
        for(int i = 0; i < 3; i++)
            R0[i] = sqrt(t0[i]*t0[i] + w0*w0);

        //line 2, pp. 1450
        Eigen::Vector3d Rplus, Rminus;
        for(int i = 0; i < 3; i++)
        {
            Rplus[i] = sqrt(tplus[i]*tplus[i] + w0*w0);
            Rminus[i] = sqrt(tminus[i]*tminus[i] + w0*w0);
        }

        // eq (11)
        Eigen::Vector3d f2;
        for(int i = 0; i < 3; i++)
        {
            double temp;
            if(w0 == 0)
            {
                if(fabs(t0[i]) < threshold)
                    temp = fabs(log(splus[i]) / sminus[i]);
                else
                    temp = (tplus[i]+splus[i]) / (tminus[i]+sminus[i]);
            }
            else
            {
                 temp = (Rplus[i]+splus[i]) / (Rminus[i]+sminus[i]);
            }
            f2[i] = log(temp);
            //fix value for points on the triangle corners
            if(f2[i] != f2[i])  //nan
                f2[i] = 0;
        }


        // eq (13) and eq (14)
        Eigen::Vector3d beta;
        double betaSum;
        if(w0 == 0)
        {
            for(int i = 0; i < 3; i++)
            {
                if(fabs(t0[i]) < threshold)
                    beta[i] = 0;
                else
                    beta[i] = atan(splus[i] / t0[i]) - atan(sminus[i] / t0[i]);
            }
        }
        else
        {
            for(int i = 0; i < 3; i++)
                beta[i] = atan((t0[i]*splus[i]) / (R0[i]*R0[i] + Rplus[i]*fabs(w0))) - atan((t0[i]*sminus[i]) / (R0[i]*R0[i] + Rminus[i]*fabs(w0)));
        }
        betaSum = beta[0] + beta[1] + beta[2];


        // eq (19), integral of kernel 1/R
        I1 = 0;
        for(int i = 0; i < 3; i++)
            I1 += t0[i]*f2[i];
        I1 -= fabs(w0) * betaSum;
    }

    __device__ static double kernel(const Eigen::Vector3d& x, const GPU_SamplePoint& sample)
    {
        double I1;
        Eigen::Vector3d Igrad;
        integrateOverTriangle(x, sample, I1, Igrad);
        //double result = (sample.value * dot(Igrad, sample.normal) - sample.derivative * I1) / (-4 * PI);
        double result = I1 * sample.derivative;
        return result;
    }

    __global__ static void kernel_integrate(const GPU_SamplePoint* samples, int* sampleCount, Eigen::Vector3d* x,
                                            double* result)
    {
        __shared__ double cache[THREADS_PER_BLOCK];
        float temp = 0;
        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        int cacheIndex = threadIdx.x;
        while(tid < *sampleCount) {
            temp += kernel(*x, samples[tid]);
            tid += blockDim.x * gridDim.x;
        }
        cache[cacheIndex] = temp;
        __syncthreads();

        int i = blockDim.x / 2;
        while(i > 0) {
            if(cacheIndex < i)
                cache[cacheIndex] += cache[cacheIndex + i];
            i /= 2;
            __syncthreads();
        }
        if(cacheIndex == 0)
            result[blockIdx.x] = cache[0];
    }

    static double CPU_getFieldValue(const GPU_SamplePoint* samples, int* sampleCount,
                                    double* gpu_result, const Eigen::Vector3d& x) {
        double cpu_result[BLOCK_COUNT];

        Eigen::Vector3d* gpu_x = nullptr;
        cudaMalloc(&gpu_x, sizeof(Eigen::Vector3d));


        cudaMemcpy(gpu_x, &x, sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);

        kernel_integrate <<<BLOCK_COUNT, THREADS_PER_BLOCK>>>(samples, sampleCount, gpu_x, gpu_result);

        cudaMemcpy(cpu_result, gpu_result, sizeof(double) * BLOCK_COUNT, cudaMemcpyDeviceToHost);

        cudaFree(gpu_x);

        double result = 0;
        for(int i = 0; i < BLOCK_COUNT; i++)
            result += cpu_result[i];
        return result;
    }

    template<class T1, class T2>
    static void toEigen(const T1& zjumat, T2& out)
    {
        for(int i = 0; i < zjumat.size(1); i++)
            for(int j = 0; j < zjumat.size(2); j++)
                out(i, j) = zjumat(i, j);
    }

    static void buildGPUSamples(const std::vector<SamplePoint>& samples, CudaVector<GPU_SamplePoint>& out)
    {
        std::vector<GPU_SamplePoint> vec;
        vec.reserve(samples.size());
        for(int i = 0; i < samples.size(); i++)
        {
            GPU_SamplePoint g;
            g.derivative = samples[i].derivative;
            toEigen(samples[i].invTransform, g.invTransform);
            toEigen(samples[i].normal, g.normal);
            toEigen(samples[i].position, g.position);
            toEigen(samples[i].transform, g.transform);
            toEigen(samples[i].tri, g.tri);
            g.size = samples[i].size;
            g.value = samples[i].value;
            vec.push_back(g);
        }
        out = CudaVector<GPU_SamplePoint>(vec);
    }

   class FieldComputerImpl
   {
   public:
       FieldComputerImpl(const std::vector<SamplePoint>& samples)
       {
           this->gpu_sampleCount = CudaPointer<int>(samples.size());
           buildGPUSamples(samples, this->gpu_samples);

           cudaMalloc(&this->gpu_integrate_result, sizeof(double) * BLOCK_COUNT);
       }

       ~FieldComputerImpl()
       {
            cudaFree(this->gpu_integrate_result);
       }

       double getFieldValue(const vec3_t &x)
       {
           Eigen::Vector3d xx;
           xx[0] = x[0];
           xx[1] = x[1];
           xx[2] = x[2];
           double result = CPU_getFieldValue(&this->gpu_samples[0], this->gpu_sampleCount.get(),
                   this->gpu_integrate_result,
                   xx);
           return result;
       }

   private:
       CudaPointer<int> gpu_sampleCount;
       CudaVector<GPU_SamplePoint> gpu_samples;
       double* gpu_integrate_result;
   };


    void FieldComputer::init(const std::vector<SamplePoint> &samples)
    {
        this->impl = new FieldComputerImpl(samples);
    }

    FieldComputer::~FieldComputer()
    {
        if(this->impl)
            delete this->impl;
    }

    double FieldComputer::getFieldValue(const vec3_t &x)
    {
        return this->impl->getFieldValue(x);
    }
}

