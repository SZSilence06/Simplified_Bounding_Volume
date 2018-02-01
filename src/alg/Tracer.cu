#include "Tracer.h"
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

    double* gpu_integrate_result;

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

        // normals of the triangle edges, Fig. 1(b)
        /*Eigen::Vector3d m[3];
        m[0][0] = v3;
        m[0][1] = l3 - u3;
        m[0][2] = 0;
        m[1][0] = -v3;
        m[1][1] = u3;
        m[1][2] = 0;
        m[2][0] = 0;
        m[2][1] = -l3;
        m[2][2] = 0;
        for(int i = 0; i < 3; i++)
            m[i] /= m[i].norm();

        // eq (34), integral of kernel grad(1/R)
        Igrad = Eigen::Vector3d::Zero();
        for(int i = 0; i < 3; i++)
            Igrad -= m[i] * f2[i];
        Eigen::Vector3d w = Eigen::Vector3d::Zero();
        w[2] = 1;
        if(w0 >= 0)
            Igrad -= betaSum * w;
        else if(w0 < 0)
            Igrad += betaSum * w;

        //transform back to world space
        Eigen::Vector4d IgradTemp;
        IgradTemp.block<3, 1>(0, 0) = Igrad;
        IgradTemp[3] = 0;
        Eigen::Vector4d IgradGlob = point.invTransform * IgradTemp;
        Igrad = IgradGlob.block<3, 1>(0, 0);*/
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

    __device__ static double getFieldValue(const GPU_SamplePoint* samples, int* sampleCount, const Eigen::Vector3d &x)
    {
        double result = 0;
        for(int i = 0; i < *sampleCount; i++)
        {
            result += kernel(x, samples[i]);
        }
        return result;
    }

    __device__ Eigen::Vector3d getGradient(const GPU_SamplePoint* samples, int* sampleCount, const Eigen::Vector3d &x)
    {
        const double STEP = 0.0001;
        double a = getFieldValue(samples, sampleCount, x);

        Eigen::Vector3d xx = x;
        xx[0] += STEP;
        double b = getFieldValue(samples, sampleCount, xx);

        Eigen::Vector3d xy = x;
        xy[1] += STEP;
        double c = getFieldValue(samples, sampleCount, xy);

        Eigen::Vector3d xz = x;
        xz[2] += STEP;
        double d = getFieldValue(samples, sampleCount, xz);

        Eigen::Vector3d result;
        result[0] = (b - a) / STEP;
        result[1] = (c - a) / STEP;
        result[2] = (d - a) / STEP;
        return result;
    }

    __global__ static void kernel_trace(const GPU_SamplePoint* samples, int* sampleCount, Eigen::Vector3d* points,
                                        Eigen::Vector3d* normals, int* pointCount, double* distance)
    {
        const double STEP = 0.01;        
        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        while(tid < *pointCount) {
            Eigen::Vector3d result = points[tid];
            Eigen::Vector3d n = normals[tid];
            double t = STEP;
            result += STEP * n;
            while(t < *distance)
            {
                Eigen::Vector3d grad = getGradient(samples, sampleCount, result);
                double norm_grad = grad.norm();
                grad /= norm_grad;
                result -= STEP * grad;
                t += STEP;
            }
            points[tid] = result;
            tid += blockDim.x * gridDim.x;
        }
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

    static double CPU_getFieldValue(const GPU_SamplePoint* samples, int* sampleCount, const Eigen::Vector3d& x) {
        double cpu_result[BLOCK_COUNT];

        Eigen::Vector3d* gpu_x = nullptr;
        cudaMalloc(&gpu_x, sizeof(Eigen::Vector3d));

        //double* gpu_result;
        //cudaMalloc(&gpu_result, sizeof(double) * BLOCK_COUNT);

        cudaMemcpy(gpu_x, &x, sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);

        kernel_integrate <<<BLOCK_COUNT, THREADS_PER_BLOCK>>>(samples, sampleCount, gpu_x, gpu_integrate_result);

        cudaMemcpy(cpu_result, gpu_integrate_result, sizeof(double) * BLOCK_COUNT, cudaMemcpyDeviceToHost);

        cudaFree(gpu_x);
        //cudaFree(gpu_result);

        double result = 0;
        for(int i = 0; i < BLOCK_COUNT; i++)
            result += cpu_result[i];
        return result;
    }

    static Eigen::Vector3d CPU_getGradient(const GPU_SamplePoint* samples, int* sampleCount, const Eigen::Vector3d &x)
    {
        const double STEP = 0.0001;
        double a = CPU_getFieldValue(samples, sampleCount, x);

        Eigen::Vector3d xx = x;
        xx[0] += STEP;
        double b = CPU_getFieldValue(samples, sampleCount, xx);

        Eigen::Vector3d xy = x;
        xy[1] += STEP;
        double c = CPU_getFieldValue(samples, sampleCount, xy);

        Eigen::Vector3d xz = x;
        xz[2] += STEP;
        double d = CPU_getFieldValue(samples, sampleCount, xz);

        Eigen::Vector3d result;
        result[0] = (b - a) / STEP;
        result[1] = (c - a) / STEP;
        result[2] = (d - a) / STEP;
        return result;
    }

    static void buildGPUVector(const matrixr_t& mat, CudaVector<Eigen::Vector3d>& out)
    {
        std::vector<Eigen::Vector3d> vec;
        vec.reserve(mat.size(2));
        for(int i = 0; i < mat.size(2); i++) {
            Eigen::Vector3d v;
            v[0] = mat(0, i);
            v[1] = mat(1, i);
            v[2] = mat(2, i);
            vec.push_back(v);
        }

        out = CudaVector<Eigen::Vector3d>(vec);
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

    /*void Tracer::tracePoints(const matrixr_t &points, const matrixr_t &normals, double length, matrixr_t &traceResult)
    {
        CudaPointer<int> gpu_sampleCount(mSamples.size());
        CudaPointer<int> gpu_pointCount(points.size(2));
        CudaPointer<double> gpu_distance(length);
        CudaVector<GPU_SamplePoint> gpu_samples;
        buildGPUSamples(mSamples, gpu_samples);

        //pass gpu_points from host to cuda kernel via std::vector
        CudaVector<Eigen::Vector3d> gpu_points, gpu_normals;
        buildGPUVector(points, gpu_points);
        buildGPUVector(normals, gpu_normals);

        kernel_trace  <<< 64, 64>>> (&gpu_samples[0], gpu_sampleCount.get(), &gpu_points[0],
                &gpu_normals[0], gpu_pointCount.get(), gpu_distance.get());
        cudaDeviceSynchronize();

        traceResult.resize(3, points.size(2));
        for(int i = 0; i < points.size(2); i++) {
            traceResult(0, i) = gpu_points[i][0];
            traceResult(1, i) = gpu_points[i][1];
            traceResult(2, i) = gpu_points[i][2];
        }
    }*/

    static void CPU_trace(const GPU_SamplePoint* samples, int* sampleCount, Eigen::Vector3d* points,
                                        Eigen::Vector3d* normals, int pointCount, double distance)
    {
        const double STEP = distance / 10;
//#pragma omp parallel for
        for(int i = 0; i < pointCount; i++) {
            std::cout << "Tracing point " << i << " / " << pointCount << std::endl;
            Eigen::Vector3d result = points[i];
            Eigen::Vector3d n = normals[i];
            double t = STEP;
            result += STEP * n;
            while(t < distance)
            {
                Eigen::Vector3d grad = CPU_getGradient(samples, sampleCount, result);
                double norm_grad = grad.norm();
                grad /= norm_grad;
                result -= STEP * grad;
                t += STEP;
            }
            points[i] = result;
        }
    }

    static void buildCPUVector(const matrixr_t& mat, std::vector<Eigen::Vector3d>& out)
    {
        out.reserve(mat.size(2));
        for(int i = 0; i < mat.size(2); i++) {
            Eigen::Vector3d v;
            v[0] = mat(0, i);
            v[1] = mat(1, i);
            v[2] = mat(2, i);
            out.push_back(v);
        }
    }

    void Tracer::tracePoints(const matrixr_t &points, const matrixr_t &normals, double length, matrixr_t &traceResult)
    {
        CudaPointer<int> gpu_sampleCount(mSamples.size());
        CudaPointer<int> gpu_pointCount(points.size(2));
        CudaPointer<double> gpu_distance(length);
        CudaVector<GPU_SamplePoint> gpu_samples;
        buildGPUSamples(mSamples, gpu_samples);

        cudaMalloc(&gpu_integrate_result, sizeof(double) * BLOCK_COUNT);

        //pass gpu_points from host to cuda kernel via std::vector
        std::vector<Eigen::Vector3d> cpu_points, cpu_normals;
        buildCPUVector(points, cpu_points);
        buildCPUVector(normals, cpu_normals);

        CPU_trace(&gpu_samples[0], gpu_sampleCount.get(), &cpu_points[0],
                &cpu_normals[0], *gpu_pointCount, *gpu_distance);

        traceResult.resize(3, points.size(2));
        for(int i = 0; i < points.size(2); i++) {
            traceResult(0, i) = cpu_points[i][0];
            traceResult(1, i) = cpu_points[i][1];
            traceResult(2, i) = cpu_points[i][2];
        }

        cudaFree(gpu_integrate_result);
    }
}

