#include "KernelRegion.h"
#include <wkylib/geometry.h>

using namespace zjucad::matrix;

namespace SBV
{
    KernelRegion::KernelRegion(const matrixr_t &points, const matrixs_t &lines, const matrixr_t& innerShell,
                               const matrixr_t& outerShell, const std::set<size_t>& innerSample,
                               const std::set<size_t>& outerSample, const TriangulatedShell& triangulation,
                               PointType collapsedPointType)
        : mPoints(points),
          mLines(lines),
          mInnerShell(innerShell),
          mOuterShell(outerShell),
          mInnerSamples(innerSample),
          mOuterSamples(outerSample),
          mTriangulation(triangulation),
          mPointType(collapsedPointType)
    {
        buildAdjacency();
        buildPolygonSequence();

        mClockwise = isClockwise();

        construct();
    }

    void KernelRegion::buildAdjacency()
    {
        //build adjacency list of the vertices.
        for(int i = 0; i < mLines.size(2); i++)
        {
            size_t a = mLines(0, i);
            size_t b = mLines(1, i);

            auto iterA = mAdjacency.find(a);
            auto iterB = mAdjacency.find(b);
            if(iterA == mAdjacency.end())
            {
                std::vector<size_t> adjacency;
                adjacency.push_back(b);
                mAdjacency.insert(std::make_pair(a, adjacency));
            }
            else
            {
                iterA->second.push_back(b);
            }

            if(iterB == mAdjacency.end())
            {
                std::vector<size_t> adjacency;
                adjacency.push_back(a);
                mAdjacency.insert(std::make_pair(b, adjacency));
            }
            else
            {
                iterB->second.push_back(a);
            }
        }
    }

    void KernelRegion::buildPolygonSequence()
    {
        //transform the polygon vertices into a sequence according to their order in the polygion.
        size_t start = mAdjacency.begin()->first;
        size_t prev = start;
        size_t next = mAdjacency.find(start)->second[0];

        mPolygon.push_back(start);
        mPolygon.push_back(next);
        while(start != next)
        {
            const std::vector<size_t>& adjacency = mAdjacency.find(next)->second;
            for(const size_t vert : adjacency)
            {
                if(vert != prev)
                {
                    prev = next;
                    next = vert;
                    mPolygon.push_back(next);
                    break;
                }
            }
        }
    }

    bool KernelRegion::isClockwise()
    {
        //judging whether the polygon is clockwise or counter-clockwise.
        double sum = 0;

        for(int i = 0; i < mPolygon.size() - 1; i++)
        {
            const matrixr_t &a = mPoints(colon(), mPolygon[i]);
            const matrixr_t &b = mPoints(colon(), mPolygon[i + 1]);

            sum += (a[1] - a[0]) * (b[1] + b[0]);
        }

        return sum > 0;
    }

    void KernelRegion::construct()
    {
        A.resize(mPolygon.size(), 3);
        for(int i = 0; i < mPolygon.size() - 1; i++)
        {
            const matrixr_t &a = mPoints(colon(), mPolygon[i]);
            const matrixr_t &b = mPoints(colon(), mPolygon[i + 1]);
            matrixr_t ab = b - a;
            ab /= norm(ab);

            if(mClockwise)
            {
                A(i, 0) = ab[1];
                A(i, 1) = -ab[0];
                A(i, 2) = ab[0] * a[1] - ab[1] * a[0];
            }
            else
            {
                A(i, 0) = -ab[1];
                A(i, 1) = ab[0];
                A(i, 2) = ab[1] * a[0] - ab[0] * a[1];
            }
        }
    }

    bool KernelRegion::contains(const matrixr_t &point)
    {
        matrixr_t homo(3, 1);
        homo[0] = point[0];
        homo[1] = point[1];
        homo[2] = 1;
        matrixr_t result = A * homo;
        for(int i = 0; i < result.size(); i++)
        {
            if(result[i] > 0)
            {
                return false;
            }
        }
        if(isInvalidRegion(point))
        {
            return false;
        }
        return true;
    }

    bool KernelRegion::isInvalidRegion(const matrixr_t &point)
    {
        for(int i = 0; i < mLines.size(2); i++)
        {
            matrixr_t triangle(2, 3);
            triangle(colon(), 0) = mPoints(colon(), mLines(0, i));
            triangle(colon(), 1) = mPoints(colon(), mLines(1, i));
            triangle(colon(), 2) = point;

            for(size_t sample : mInnerSamples)
            {
                matrixr_t bary;
                if(WKYLIB::barycentric_2D(mInnerShell(colon(), sample), triangle, bary))
                {
                    //the point is inside the tetrahedron
                    double f0 = mTriangulation.getFValue(mLines(0, i));
                    double f1 = mTriangulation.getFValue(mLines(1, i));
                    double f2 = mTriangulation.getFValue(mPointType);

                    double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2];
                    if(f > 0)
                    {
                        return true;
                    }
                }
            }

            for(size_t sample : mOuterSamples)
            {
                matrixr_t bary;
                if(WKYLIB::barycentric_2D(mOuterShell(colon(), sample), triangle, bary))
                {
                    //the point is inside the tetrahedron
                    double f0 = mTriangulation.getFValue(mLines(0, i));
                    double f1 = mTriangulation.getFValue(mLines(1, i));
                    double f2 = mTriangulation.getFValue(mPointType);

                    double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2];
                    if(f < 0)
                    {
                        return true;
                    }
                }
            }
        }

        return false;
    }
}