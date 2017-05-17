#ifndef WKY_KERNEL_REGION_H
#define WKY_KERNEL_REGION_H

#include "Common.h"
#include <map>
#include <set>
#include "TriangulatedShell.h"
#include "Shell.h"

namespace SBV
{
    class KernelRegion
    {
    public:
        KernelRegion(const std::vector<Eigen::Vector2i>& lines, const Shell& shell,
                     const std::set<size_t>& innerSample, const std::set<size_t>& outerSample,
                     const TriangulatedShell& triangulation, PointType collapsedPointType);

        bool contains(const Point& point) const;
        bool isInvalidRegion(const Point& point) const;

    private:                      
        void buildAdjacency();
        void buildPolygonSequence();
        bool isClockwise();
        void construct();

    private:
        const std::vector<Eigen::Vector2i>& mLines;
        const Shell& mShell;
        const std::set<size_t>& mInnerSamples;
        const std::set<size_t>& mOuterSamples;
        const TriangulatedShell& mTriangulation;
        PointType mPointType;
        bool mClockwise;

        std::map<size_t, std::vector<size_t> > mAdjacency;
        std::vector<size_t> mPolygon;   //recording the polygon verts, in a cycled sequence.

        Eigen::MatrixXd A;

        friend class SamplingQuadTree;
    };
}

#endif
