#ifndef WKY_SAMPLING_QUAD_TREE_H
#define WKY_SAMPLING_QUAD_TREE_H

#include "Common.h"
#include <wkylib/Geometry/QuadTree.h>

namespace SBV
{
    class KernelRegion;

    class SamplingQuadTree : private WKYLIB::Geometry::QuadTree
    {
    public:
        SamplingQuadTree(const KernelRegion& kernel, double xmax, double xmin, double ymax, double ymin, double sampleRadius);

    private:
        void build();
        void sample(QuadTreeNode* node);
        bool isOutsideKernelRegion(QuadTreeNode* node);

    private:
        const KernelRegion& mKernel;
        double mSampleRadius;
    };
}

#endif
