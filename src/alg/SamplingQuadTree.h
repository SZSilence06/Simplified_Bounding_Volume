#ifndef WKY_SAMPLING_QUAD_TREE_H
#define WKY_SAMPLING_QUAD_TREE_H

#include "Common.h"

namespace SBV
{
    class KernelRegion;

    /**
     * @brief The class to sample a region and keep those inside a kernel region.
     *
     * This class is so far not used. It does not use a quad tree to sample. Just performs the simple uniform sample strategy so far.
     */
    class SamplingQuadTree
    {
    public:
        SamplingQuadTree(const KernelRegion& kernel, double xmax, double xmin, double ymax, double ymin,
                         double zmax, double zmin, double sampleRadius);

        inline const std::vector<matrixr_t>& getSamples() const { return mSamples; }

    private:
        void build();
        void sample(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
        bool isOutsideKernelRegion(double xmin, double xmax, double ymin, double ymax);

    private:
        const KernelRegion& mKernel;
        double mSampleRadius;
        std::vector<matrixr_t> mSamples;
    };
}

#endif
