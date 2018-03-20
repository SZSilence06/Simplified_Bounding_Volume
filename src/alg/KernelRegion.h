#ifndef WKY_KERNEL_REGION_H
#define WKY_KERNEL_REGION_H

#include "Common.h"
#include <map>
#include <set>
#include "TriangulatedShell.h"
#include "Shell.h"

namespace SBV
{
     /**
     * @brief The KernelRegion class
     *
     * This class describes the kernel region using a set of linear equations.
     * These linear equations take the following form : Ax < b, where A is a matrix and x,b are vectors.
     */
    class KernelRegion
    {
    public:
        /**
         * @brief KernelRegion constructor
         * @param points : point matrix of the boundary of the one ring area
         * @param faces : triangle matrix of the boundary of the one ring area
         * @param related_vert_for_boundary_faces : the corresponding vertex id on the collapsed edge for each triangle on the boundary of the one ring area.
         *                                          This parameter is used for determing the sign of each linear equation.
         * @param shell : the input samples
         * @param innerSample : the id of samples that are on outer shell inside the kernel region.
         * @param outerSample : the id of samples that are on outer shell inside the kernel region.
         * @param triangulation : the input delaunay triangulation
         * @param collapsedPointType : point type of the collapsed edge(inner shell point, outer shell point or zero set point).
         */
        KernelRegion(const matrixr_t& points, const matrixs_t& faces,
                     const std::vector<size_t>& related_vert_for_boundary_faces,
                     const Shell& shell, const std::set<size_t>& innerSample, const std::set<size_t>& outerSample,
                     const TriangulatedShell& triangulation, PointType collapsedPointType);

        /**
         * @brief test whether a given point is inside the kernel region.
         * @param point : the point to test.
         * @return true if inside and false otherwise.
         */
        bool contains(const vec3_t& point) const;

    private:                      
        void construct(const std::vector<size_t>& related_vert_for_boundary_faces);

        bool isInvalidRegion(const vec3_t& point) const;

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
