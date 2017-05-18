#ifndef WKY_SAMPLING_QUAD_TREE_H
#define WKY_SAMPLING_QUAD_TREE_H

#include "Common.h"

namespace SBV
{
    class KernelRegion;

    class SamplingQuadTree
    {
    public:
#ifdef VER_2D
        SamplingQuadTree(const KernelRegion& kernel, double xmax, double xmin, double ymax, double ymin, double sampleRadius);
#else
        SamplingQuadTree(const KernelRegion& kernel, double xmax, double xmin, double ymax, double ymin,
                         double zmax, double zmin, double sampleRadius);
#endif

        inline const std::vector<Point>& getSamples() const { return mSamples; }

    private:
#ifdef VER_2D
        void sample(double xmin, double xmax, double ymin, double ymax);
#else
        void sample(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
#endif

    private:
        const KernelRegion& mKernel;
        double mSampleRadius;
        std::vector<Point> mSamples;
    };
}

#endif
