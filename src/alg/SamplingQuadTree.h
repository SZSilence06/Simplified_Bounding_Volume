#ifndef WKY_SAMPLING_QUAD_TREE_H
#define WKY_SAMPLING_QUAD_TREE_H

#include "Common.h"

namespace SBV
{
    class KernelRegion;

    class SamplingQuadTree
    {
    public:
        SamplingQuadTree(const KernelRegion& kernel, double xmax, double xmin, double ymax, double ymin, double sampleRadius);

        inline const std::vector<matrixr_t>& getSamples() const { return mSamples; }

    private:
        void build();
        void sample(double xmin, double xmax, double ymin, double ymax);
        bool isOutsideKernelRegion(double xmin, double xmax, double ymin, double ymax);

    private:
        const KernelRegion& mKernel;
        double mSampleRadius;
        std::vector<matrixr_t> mSamples;
    };
}

#endif
