
#include "SamplingQuadTree.h"
#include "KernelRegion.h"
#include "eigen3.3/Eigen/Dense"
#include <iostream>

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
        int x = 0;
        int y = 0;
        int xCount = (xmax - xmin) / mSampleRadius;
        int yCount = (ymax - ymin) / mSampleRadius;

        while(x < xCount && y < yCount)
        {
            Point point;
            point[0] = xmin + mSampleRadius * x;
            point[1] = ymin + mSampleRadius * y;

            if(mKernel.contains(point))
            {
                mSamples.push_back(point);
            }
            x++;
            if(x >= xCount)
            {
                x = 0;
                y++;
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
