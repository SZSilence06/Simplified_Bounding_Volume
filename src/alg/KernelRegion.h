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
        KernelRegion(const matrixr_t& points, const matrixs_t& faces,
                     const std::vector<size_t>& related_vert_for_boundary_faces,
                     const Shell& shell, const std::set<size_t>& innerSample, const std::set<size_t>& outerSample,
                     const TriangulatedShell& triangulation, PointType collapsedPointType);

        bool contains(const vec3_t& point) const;
        bool isInvalidRegion(const vec3_t& point) const;

    private:                      
        void construct(const std::vector<size_t>& related_vert_for_boundary_faces);

    private:
        const matrixr_t& mPoints;
        const matrixs_t& mFaces;
        const Shell& mShell;
        const std::set<size_t>& mInnerSamples;
        const std::set<size_t>& mOuterSamples;
        const TriangulatedShell& mTriangulation;
        PointType mPointType;
        bool mClockwise;

        std::map<size_t, std::vector<size_t> > mAdjacency;
        std::vector<size_t> mPolygon;   //recording the polygon verts, in a cycled sequence.

        matrixr_t A;

        friend class SamplingQuadTree;
    };
}

#endif
