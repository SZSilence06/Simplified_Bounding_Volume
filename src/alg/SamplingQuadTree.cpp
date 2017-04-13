
#include "SamplingQuadTree.h"
#include "KernelRegion.h"

namespace SBV
{
    SamplingQuadTree::SamplingQuadTree(const KernelRegion &kernel, double xmax, double xmin, double ymax, double ymin,
                                       double sampleRadius)
        : QuadTree(xmin, xmax, ymin, ymax),
          mKernel(kernel),
          mSampleRadius(sampleRadius)
    {
        build();
    }

    void SamplingQuadTree::build()
    {
        sample(this->mRoot);
    }

    void SamplingQuadTree::sample(QuadTreeNode *node)
    {
        matrixr_t lb(2, 1);
        matrixr_t lt(2, 1);
        matrixr_t rb(2, 1);
        matrixr_t rt(2, 1);

        lb[0] = node->xMin;
        lb[1] = node->yMin;
        lt[0] = node->xMin;
        lt[1] = node->yMax;
        rb[0] = node->xMax;
        rb[1] = node->yMin;
        rt[0] = node->xMax;
        rt[1] = node->yMax;

        bool invalidLB = mKernel.isInvalidRegion(lb);
        bool invalidLT = mKernel.isInvalidRegion(lt);
        bool invalidRB = mKernel.isInvalidRegion(rb);
        bool invalidRT = mKernel.isInvalidRegion(rt);

        if(isOutsideKernelRegion(node))
        {
            return;
        }

        if(node->xMax - node->xMin < mSampleRadius || node->yMax - node->yMin < mSampleRadius)
        {
            matrixr_t point(2, 1);
            point[0] = (node->xMax + node->xMin) / 2;
            point[1] = (node->yMax + node->yMin) / 2;
            insertToNode(point, node);
            return;
        }

        split(node);
        sample(node->lb);
        sample(node->lt);
        sample(node->rb);
        sample(node->rt);
    }
}
