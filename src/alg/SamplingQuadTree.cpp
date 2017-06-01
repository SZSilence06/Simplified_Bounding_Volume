
#include "SamplingQuadTree.h"
#include "KernelRegion.h"
#include <iostream>
#include <mutex>

namespace SBV
{
#ifdef VER_2D
    SamplingQuadTree::SamplingQuadTree(const KernelRegion &kernel, double xmax, double xmin, double ymax, double ymin,
                                       double sampleRadius)
        : mKernel(kernel),
          mSampleRadius(sampleRadius)
    {
        sample(xmin, xmax, ymin, ymax);
    }
#else
    SamplingQuadTree::SamplingQuadTree(const KernelRegion &kernel, double xmax, double xmin, double ymax, double ymin,
                                       double zmax, double zmin, double sampleRadius)
        : mKernel(kernel),
          mSampleRadius(sampleRadius)
    {
        sample(xmin, xmax, ymin, ymax, zmax, zmin);
    }
#endif

#ifdef VER_2D
    void SamplingQuadTree::sample(double xmin, double xmax, double ymin, double ymax)
    {
        const int xCount = (xmax - xmin) / mSampleRadius;
        const int yCount = (ymax - ymin) / mSampleRadius;

        std::mutex mtx;

#pragma omp parallel for schedule(dynamic, 1)
        for(int x = 0; x < xCount; x++)
        {
            for(int y = 0; y < yCount; y++)
            {
                Point point;
                point[0] = xmin + mSampleRadius * x;
                point[1] = ymin + mSampleRadius * y;

                if(mKernel.contains(point))
                {
                    mtx.lock();
                    mSamples.push_back(point);
                    mtx.unlock();
                }
            }
        }
    }
#else
    void SamplingQuadTree::sample(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
    {
        int x = 0;
        int y = 0;
        int z = 0;
        int xCount = (xmax - xmin) / mSampleRadius;
        int yCount = (ymax - ymin) / mSampleRadius;
        int zCount = (zmax - zmin) / mSampleRadius;

        while(x < xCount && y < yCount && z < zCount)
        {
            Point point;
            point[0] = xmin + mSampleRadius * x;
            point[1] = ymin + mSampleRadius * y;
            point[2] = zmin + mSampleRadius * z;

            if(mKernel.contains(point))
            {
                mSamples.push_back(point);
            }
            x++;
            if(x >= xCount)
            {
                x = 0;
                y++;
                if(y >= yCount)
                {
                    y = 0;
                    z++;
                }
            }
        }
    }
#endif
}
