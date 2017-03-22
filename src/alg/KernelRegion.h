#ifndef WKY_KERNEL_REGION_H
#define WKY_KERNEL_REGION_H

#include "Common.h"
#include <map>
#include <set>

namespace SBV
{
    class KernelRegion
    {
    public:
        KernelRegion(const matrixr_t& points, const matrixs_t& lines, const matrixr_t& innerShell, const matrixr_t& outerShell,
                     const std::set<size_t>& innerSample, const std::set<size_t>& outerSample);

        bool contains(const matrixr_t& point);

    private:                      
        void buildAdjacency();
        void buildPolygonSequence();
        bool isClockwise();
        void construct();

    private:
        const matrixr_t& mPoints;
        const matrixs_t& mLines;
        const matrixr_t& mInnerShell;
        const matrixr_t& mOuterShell;
        const std::set<size_t>& mInnerSamples;
        const std::set<size_t>& mOuterSamples;
        bool mClockwise;

        std::map<size_t, std::vector<size_t> > mAdjacency;
        std::vector<size_t> mPolygon;   //recording the polygon verts, in a cycled sequence.

        matrixr_t A;
    };
}

#endif
