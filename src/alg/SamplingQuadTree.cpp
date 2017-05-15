
#include "SamplingQuadTree.h"
#include "KernelRegion.h"
#include <iostream>

namespace SBV
{
    SamplingQuadTree::SamplingQuadTree(const KernelRegion &kernel, double xmax, double xmin, double ymax, double ymin,
                                       double sampleRadius)
        : mKernel(kernel),
          mSampleRadius(sampleRadius)
    {
        sample(xmin, xmax, ymin, ymax);
    }

    void SamplingQuadTree::sample(double xmin, double xmax, double ymin, double ymax)
    {
        int x = 0;
        int y = 0;
        int xCount = (xmax - xmin) / mSampleRadius;
        int yCount = (ymax - ymin) / mSampleRadius;

        while(x < xCount && y < yCount)
        {
            matrixr_t point(2, 1);
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
        /*if(isOutsideKernelRegion(xmin, xmax, ymin, ymax))
        {
            //return;
        }

        double xmid = (xmin + xmax) / 2;
        double ymid = (ymin + ymax) / 2;

        if(xmax - xmin < mSampleRadius || ymax - ymin < mSampleRadius)
        {
            matrixr_t point(2, 1);
            point[0] = xmid;
            point[1] = ymid;
            if(mKernel.contains(point))
            {
                mSamples.push_back(point);
            }
            return;
        }
        sample(xmin, xmid, ymin, ymid);
        sample(xmin, xmid, ymid, ymax);
        sample(xmid, xmax, ymin, ymid);
        sample(xmid, xmax, ymid, ymax);*/
    }

    bool SamplingQuadTree::isOutsideKernelRegion(double xmin, double xmax, double ymin, double ymax)
    {
        matrixr_t p1(2, 1), p2(2, 1), p3(2, 1), p4(2, 1);
        p1[0] = xmin;
        p1[1] = ymin;
        p2[0] = xmin;
        p2[1] = ymax;
        p3[0] = xmax;
        p3[1] = ymin;
        p4[0] = xmax;
        p4[1] = ymax;

        return false;

        if(mKernel.contains(p1) || mKernel.contains(p2) || mKernel.contains(p3) || mKernel.contains(p4))
        {
            return false;
        }
    }
}
