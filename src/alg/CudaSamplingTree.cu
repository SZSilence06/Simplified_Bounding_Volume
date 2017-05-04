#include "CudaSamplingTree.h"
#include "CudaKernelRegion.h"

using namespace WKYLIB::Cuda;

namespace SBV
{
    __global__ void computeSample(CudaPointer<CudaKernelRegion> kernel,
                                  CudaPointer<double> xmin,
                                  CudaPointer<double> ymin,
                                  CudaPointer<int> xCount,
                                  CudaPointer<int> yCount,
                                  CudaPointer<double> sampleRadius,
                                  CudaPointer<CudaVector<Eigen::Vector2d>> samples)
    {
        int x = threadIdx.x + blockIdx.x * blockDim.x;
        int y = threadIdx.y + blockIdx.y * blockDim.y;

        while(x < *xCount && y < *yCount)
        {
            Eigen::Vector2d point;
            point[0] = *xmin + *sampleRadius * x;
            point[1] = *ymin + *sampleRadius * y;

            if(kernel->contains(point))
            {
                samples->push_back(point);
            }

            x += blockDim.x * gridDim.x;
            y += blockDim.y * gridDim.y;
        }
    }

    CudaSamplingTree::CudaSamplingTree(CudaPointer<CudaKernelRegion> kernel, double xmax, double xmin, double ymax, double ymin,
                                       double sampleRadius)
        : mKernel(kernel),
          mSampleRadius(sampleRadius)
    {
        sample(xmin, xmax, ymin, ymax);
    }

    void CudaSamplingTree::sample(double xmin, double xmax, double ymin, double ymax)
    {
        int xCount = (xmax - xmin) / mSampleRadius;
        int yCount = (ymax - ymin) / mSampleRadius;

        CudaPointer<int> gpu_xCount(xCount);
        CudaPointer<int> gpu_yCount(yCount);
        CudaPointer<double> gpu_xmin(xmin);
        CudaPointer<double> gpu_ymin(ymin);
        CudaPointer<double> gpu_sampleRadius(mSampleRadius);

        dim3 grid(256, 256);
        computeSample <<<grid, grid>>> (mKernel, gpu_xmin, gpu_ymin, gpu_xCount,
                                     gpu_yCount, gpu_sampleRadius, mSamples);

        cudaDeviceSynchronize();
    }
}
