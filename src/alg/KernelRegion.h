#ifndef WKY_KERNEL_REGION_H
#define WKY_KERNEL_REGION_H

#include "Common.h"
#include <map>
#include <set>
#include "TriangulatedShell.h"
#include "Shell.h"


namespace SBV
{
    enum InvalidRegionType
    {
        INVALID_REGION_BOUNDARY,
        INVALID_REGION_ZERO_SET
    };

    class KernelRegion
    {
    public:
        KernelRegion(const matrixr_t& points, const matrixs_t& lines, const Shell& shell,
                     const std::set<size_t>& innerSample, const std::set<size_t>& outerSample,
                     size_t maxErrorInnerSample, size_t maxErrorOuterSample,
                     const TriangulatedShell& triangulation, PointType collapsedPointType,
                     InvalidRegionType invalidRegionType);

        bool contains(const matrixr_t& point) const;
        bool isInvalidRegion(const matrixr_t& point) const;
        bool getBestPos(const matrixr_t& Q, matrixr_t& output_position) const;

    private:                      
        void buildAdjacency();
        void buildPolygonSequence();
        bool isClockwise();
        void construct();
        void computeInvalidConstraints();
        void buildSegment(const matrixr_t& a, const matrixr_t& b, matrixr_t& segment);
        void buildConstraintForBundary(const matrixr_t& a, const matrixr_t& b, const matrixr_t& E, int situation);

    private:
        const matrixr_t& mPoints;
        const matrixs_t& mLines;
        const Shell& mShell;
        const std::set<size_t>& mInnerSamples;
        const std::set<size_t>& mOuterSamples;
        size_t mInnerMaxErrorSample;
        size_t mOuterMaxErrorSample;
        const TriangulatedShell& mTriangulation;
        PointType mPointType;
        InvalidRegionType mInvalidRegionType;
        bool mClockwise;

        std::map<size_t, std::vector<size_t> > mAdjacency;
        std::vector<size_t> mPolygon;   //recording the polygon verts, in a cycled sequence.

        matrixr_t A;
        std::vector<matrixr_t> mInvalidConstraints;

        friend class SamplingQuadTree;
    };
}

#endif
