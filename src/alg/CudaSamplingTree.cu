#include "CudaSamplingTree.h"
#include "CudaKernelRegion.h"

using namespace WKYLIB::Cuda;

namespace SBV
{
    __global__ void computeSample(CudaKernelRegion* kernel,
                                  double* xmin,
                                  double* ymin,
                                  int* xCount,
                                  int* yCount,
                                  double* sampleRadius,
                                  CudaVector<Eigen::Vector2d>* samples)
    {
        int originalX = threadIdx.x + blockIdx.x * blockDim.x;
        int x = originalX;
        int y = threadIdx.y + blockIdx.y * blockDim.y;

        //kernel->printTest();

        while(x < *xCount && y < *yCount)
        {
            Eigen::Vector2d point;
            point[0] = *xmin + *sampleRadius * x;
            point[1] = *ymin + *sampleRadius * y;        

            if(kernel->contains(point))
            {
                samples->push_back_on_device(point);
            }

            x += blockDim.x * gridDim.x;
            if(x >= *xCount)
            {
                x = originalX;
                y += blockDim.y * gridDim.y;
            }
        }
    }

    /*void computeSample_cpu(CudaKernelRegion* kernel,
                                  double* xmin,
                                  double* ymin,
                                  int* xCount,
                                  int* yCount,
                                  double* sampleRadius,
                                  CudaVector<Eigen::Vector2d>* samples)
    {
        int x = 0;
        int y = 0;

        //kernel->printTest();

        while(x < *xCount && y < *yCount)
        {
            Eigen::Vector2d point;
            point[0] = *xmin + *sampleRadius * x;
            point[1] = *ymin + *sampleRadius * y;

            if(kernel->contains(point))
            {
                samples->push_back(point);
            }

            x++;
            y++;
        }
    }*/

    CudaSamplingTree::CudaSamplingTree(CudaPointer<CudaKernelRegion> kernel, double xmax, double xmin, double ymax, double ymin,
                                       double sampleRadius)
        : mKernel(kernel),
          mSampleRadius(sampleRadius)
    {
        mSamples.assign(CudaVector<Eigen::Vector2d>());
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

        mSamples->reserve(xCount * yCount);

        computeSample <<<1, 1>>> (mKernel.get(), gpu_xmin.get(), gpu_ymin.get(), gpu_xCount.get(),
                                     gpu_yCount.get(), gpu_sampleRadius.get(), mSamples.get());

        //computeSample_cpu(mKernel.get(), gpu_xmin.get(), gpu_ymin.get(), gpu_xCount.get(),
        //                  gpu_yCount.get(), gpu_sampleRadius.get(), mSamples.get());

        cudaDeviceSynchronize();
    }
}
