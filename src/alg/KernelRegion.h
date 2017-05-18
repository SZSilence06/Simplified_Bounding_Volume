#ifndef WKY_KERNEL_REGION_H
#define WKY_KERNEL_REGION_H

#include "Common.h"
#include <map>
#include <set>
#include "TriangulatedShell.h"
#include "InvalidRegionType.h"

namespace SBV
{
    class Shell;

    class KernelRegion
    {
    public:
        KernelRegion(const std::vector<Eigen::Vector2i>& lines, const Shell& shell,
                     const std::set<size_t>& innerSample, const std::set<size_t>& outerSample,
                     const TriangulatedShell& triangulation, PointType collapsedPointType,
                     InvalidRegionType invalidRegionType);

        bool contains(const Point& point) const;
        bool isInvalidRegion(const Point& point) const;

    private:                      
        void buildAdjacency();
        void buildPolygonSequence();
        bool isClockwise();
        void construct();

        //void computeInvalidConstraints();
        //void buildSegment(const matrixr_t& a, const matrixr_t& b, matrixr_t& segment);
        //void buildConstraintForBundary(const matrixr_t& a, const matrixr_t& b, const matrixr_t& E, int situation);

    private:
        const std::vector<Eigen::Vector2i>& mLines;
        const Shell& mShell;
        const std::set<size_t>& mInnerSamples;
        const std::set<size_t>& mOuterSamples;
        const TriangulatedShell& mTriangulation;
        PointType mPointType;
        InvalidRegionType mInvalidRegionType;
        bool mClockwise;

        std::map<size_t, std::vector<size_t> > mAdjacency;
        std::vector<size_t> mPolygon;   //recording the polygon verts, in a cycled sequence.


        //std::vector<matrixr_t> mInvalidConstraints;
        Eigen::MatrixXd A;

        friend class SamplingQuadTree;
        friend class CudaKernelRegion;
    };
}

#endif
