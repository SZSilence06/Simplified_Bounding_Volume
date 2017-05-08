#include "KernelRegion.h"
#include "Shell.h"
#include <wkylib/geometry.h>

using namespace zjucad::matrix;

namespace SBV
{
    KernelRegion::KernelRegion(const matrixr_t &points, const matrixs_t &lines, const Shell& shell,
                               const std::set<size_t>& innerSample, const std::set<size_t>& outerSample,
                               const TriangulatedShell& triangulation, PointType collapsedPointType,
                               InvalidRegionType invalidRegionType)
        : mPoints(points),
          mLines(lines),
          mShell(shell),
          mInnerSamples(innerSample),
          mOuterSamples(outerSample),
          mTriangulation(triangulation),
          mPointType(collapsedPointType),
          mInvalidRegionType(invalidRegionType)
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

        //computeInvalidConstraints();
    }

    void KernelRegion::computeInvalidConstraints()
    {
        if(mInvalidRegionType == INVALID_REGION_BOUNDARY)
        {
            for(size_t sample : mInnerSamples)
            {
                 matrixr_t EInner = mShell.mInnerShell(colon(), sample);
                 for(int i = 0; i < mLines.size(2); i++)
                 {
                     const matrixr_t& a = mPoints(colon(), mLines(0, i));
                     const matrixr_t& b = mPoints(colon(), mLines(1, i));
                     double signA = mTriangulation.getSign(mLines(0, i));
                     double signB = mTriangulation.getSign(mLines(1, i));

                     if(mPointType == POINT_INNER)
                     {
                         //F(T) = F(E), so it is situation 1 or 3
                         if(signA < 0 && signB < 0)
                         {

                         }
                         else if(signA > 0 && signB > 0)
                         {
                             //situation 1
                             buildConstraintForBundary(a, b, EInner, 1);
                         }
                         else
                         {
                             //situation 3
                             if(signA > 0)
                             {
                                 buildConstraintForBundary(a, b, EInner, 3);
                             }
                             else
                             {
                                 buildConstraintForBundary(b, a, EInner, 3);
                             }
                         }
                     }
                     else
                     {
                         //F(T) != F(E), so it is situation 2 or 4
                         if(signA < 0 && signB < 0)
                         {
                             buildConstraintForBundary(a, b, EInner, 2);
                         }
                         else if(signA > 0 && signB > 0)
                         {

                         }
                         else
                         {
                             if(signA > 0)
                             {
                                 buildConstraintForBundary(b, a, EInner, 4);
                             }
                             else
                             {
                                 buildConstraintForBundary(a, b, EInner, 4);
                             }
                         }
                     }
                 }
            }

            for(size_t sample : mOuterSamples)
            {
                 matrixr_t EOuter = mShell.mOuterShell(colon(), sample);
                 for(int i = 0; i < mLines.size(2); i++)
                 {
                     const matrixr_t& a = mPoints(colon(), mLines(0, i));
                     const matrixr_t& b = mPoints(colon(), mLines(1, i));
                     double signA = mTriangulation.getSign(mLines(0, i));
                     double signB = mTriangulation.getSign(mLines(1, i));

                     if(mPointType == POINT_OUTER)
                     {
                         //F(T) > 0
                          //F(T) = F(E), so it is situation 1 or 3
                         if(signA < 0 && signB < 0)
                         {
                             buildConstraintForBundary(a, b, EOuter, 1);
                         }
                         else if(signA > 0 && signB > 0)
                         {

                         }
                         else
                         {
                             if(signA > 0)
                             {
                                 buildConstraintForBundary(b, a, EOuter, 3);
                             }
                             else
                             {
                                 buildConstraintForBundary(a, b, EOuter, 3);
                             }
                         }
                     }
                     else
                     {
                         //F(T) < 0
                          //F(T) != F(E), so it is situation 2 or 4
                         if(signA < 0 && signB < 0)
                         {

                         }
                         else if(signA > 0 && signB > 0)
                         {
                             buildConstraintForBundary(a, b, EOuter, 2);
                         }
                         else
                         {
                             if(signA > 0)
                             {
                                 buildConstraintForBundary(a, b, EOuter, 4);
                             }
                             else
                             {
                                 buildConstraintForBundary(b, a, EOuter, 4);
                             }
                         }
                     }
                 }
            }

            /*matrixr_t EInner = mShell.mInnerShell(colon(), mInnerMaxErrorSample);
            matrixr_t EOuter = mShell.mOuterShell(colon(), mOuterMaxErrorSample);
            for(int i = 0; i < mLines.size(2); i++)
            {
                const matrixr_t& a = mPoints(colon(), mLines(0, i));
                const matrixr_t& b = mPoints(colon(), mLines(1, i));

                if(mTriangulation.getSign(mLines(0, i)) < 0 && mTriangulation.getSign(mLines(1, i)) < 0)
                {
                    //the vertices of the edge has the same sign with E, so no invalid region
                    buildConstraintForBundary(a, b, EOuter, 1);
                }
                else if(mTriangulation.getSign(mLines(0, i)) > 0 && mTriangulation.getSign(mLines(1, i)) > 0)
                {
                    buildConstraintForBundary(a, b, EInner, 1);
                }
                else
                {
                    buildConstraintForBundary(a, b, EInner, 2);
                    buildConstraintForBundary(a, b, EOuter, 2);
                }
            }*/
        }
    }

    void KernelRegion::buildConstraintForBundary(const matrixr_t &a, const matrixr_t &b, const matrixr_t &E, int situation)
    {
        matrixr_t constraint(3, 3);
        matrixr_t constraint_aE, constraint_bE;

        buildSegment(a, E, constraint_aE);
        if(b[0] * constraint_aE[0] + b[1] * constraint_aE[1] + constraint_aE[2] > 0)
        {
            constraint(0, colon()) = constraint_aE;
        }
        else
        {
            constraint(0, colon()) = -constraint_aE;
        }

        buildSegment(b, E, constraint_bE);
        if(a[0] * constraint_bE[0] + a[1] * constraint_bE[1] + constraint_bE[2] > 0)
        {
            constraint(1, colon()) = constraint_bE;
        }
        else
        {
            constraint(1, colon()) = -constraint_bE;
        }

        if(situation < 3)
        {
            //situation 1, 2
            matrixr_t constraint_ab;
            buildSegment(a, b, constraint_ab);
            double dist = E[0] * constraint_ab[0] + E[1] * constraint_ab[1] + constraint_ab[2];
            constraint_ab[2] -= 2 * dist;
            if(situation == 1)
            {
                if(dist < 0)
                {
                    constraint(2, colon()) = constraint_ab;
                }
                else
                {
                    constraint(2, colon()) = -constraint_ab;
                }
            }
            else
            {
                if(dist > 0)
                {
                    constraint(2, colon()) = constraint_ab;
                }
                else
                {
                    constraint(2, colon()) = -constraint_ab;
                }
            }
        }
        else
        {
            //situation 3, 4
            matrixr_t Y = (a + b) / 2;
            matrixr_t constraint_EY;
            buildSegment(E, Y, constraint_EY);
            double dist = a[0] * constraint_EY[0] + a[1] * constraint_EY[1] + constraint_EY[2];
            constraint_EY[2] += dist;
            if(situation == 3)
            {
                if(dist > 0)
                {
                    constraint(2, colon()) = constraint_EY;
                }
                else
                {
                    constraint(2, colon()) = -constraint_EY;
                }
            }
            else
            {
                if(dist < 0)
                {
                    constraint(2, colon()) = constraint_EY;
                }
                else
                {
                    constraint(2, colon()) = -constraint_EY;
                }
            }
        }

        mInvalidConstraints.push_back(constraint);
    }

    void KernelRegion::buildSegment(const matrixr_t &a, const matrixr_t &b, matrixr_t& segment)
    {
        matrixr_t ab = b - a;
        ab /= norm(ab);
        segment.resize(1, 3);
        segment[0] = ab[1];
        segment[1] = -ab[0];
        segment[2] = ab[0] * a[1] - ab[1] * a[0];
    }

    bool KernelRegion::contains(const matrixr_t &point) const
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

    bool KernelRegion::isInvalidRegion(const matrixr_t &point) const
    {
        /*bool result1 = false;
        if(mInvalidRegionType == INVALID_REGION_BOUNDARY)
        {
            for(int i = 0; i < mInvalidConstraints.size(); i++)
            {
                bool invalid = true;
                const matrixr_t& constraint = mInvalidConstraints[i];
                matrixr_t homo(3, 1);
                homo[0] = point[0];
                homo[1] = point[1];
                homo[2] = 1;
                matrixr_t result = constraint * homo;
                for(int j = 0; j < result.size(); j++)
                {
                    if(result[j] > 0)
                    {
                        invalid = false;
                        break;
                    }
                }
                if(invalid)
                {
                    result1 = true;
                    break;
                }
            }
            //return false;
        }*/
        for(int i = 0; i < mLines.size(2); i++)
        {
            matrixr_t triangle(2, 3);
            triangle(colon(), 0) = mPoints(colon(), mLines(0, i));
            triangle(colon(), 1) = mPoints(colon(), mLines(1, i));
            triangle(colon(), 2) = point;

            for(size_t sample : mInnerSamples)
            {
                matrixr_t bary;
                if(WKYLIB::barycentric_2D(mShell.mInnerShell(colon(), sample), triangle, bary))
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
                if(WKYLIB::barycentric_2D(mShell.mOuterShell(colon(), sample), triangle, bary))
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

    bool KernelRegion::getBestPos(const matrixr_t &Q, matrixr_t &output_position) const
    {

    }
}
