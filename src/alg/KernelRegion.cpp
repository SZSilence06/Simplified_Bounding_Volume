#include "KernelRegion.h"
#include "BaryComputer.h"
#include <wkylib/geometry.h>
#include <iostream>
#include <zjucad/matrix/io.h>

using namespace zjucad::matrix;

namespace SBV
{
    KernelRegion::KernelRegion(const matrixr_t& points, const matrixs_t& faces, const matrixr_t& onePointInRegion,
                               const Shell& shell, const std::set<size_t>& innerSample, const std::set<size_t>& outerSample,
                               const TriangulatedShell& triangulation, PointType collapsedPointType)
        : mPoints(points),
          mFaces(faces),
          mShell(shell),
          mInnerSamples(innerSample),
          mOuterSamples(outerSample),
          mTriangulation(triangulation),
          mPointType(collapsedPointType)
    {
        construct(onePointInRegion);
    }

    void KernelRegion::construct(const matrixr_t& onePointInRegion)
    {
        A.resize(mFaces.size(2), 4);
        for(int i = 0; i < mFaces.size(2); i++)
        {
            const vec3_t a = mPoints(colon(), mFaces(0, i));
            const vec3_t b = mPoints(colon(), mFaces(1, i));
            const vec3_t c = mPoints(colon(), mFaces(2, i));

            vec3_t n = cross(b- a, c- a);
            n /= norm(n);
            if(dot(a - onePointInRegion, n) < 0)
            {
                n = -n;
            }

            A(i, 0) = n[0];
            A(i, 1) = n[1];
            A(i, 2) = n[2];
            A(i, 3) = -n[0]*a[0] - n[1]*a[1] - n[2]*a[2];
        }
    }

    /*void KernelRegion::findShellSamples()
    {
        mInnerSamples.clear();
        mOuterSamples.clear();

        matrixr_t polygon(2, mLines.size(2));
        for(int i = 0; i < mLines.size(2); i++)
        {
            polygon(colon(), i) = mPoints(colon(), mPolygon[i]);
        }

        matrixs_t points;
        mShell.getInnerTree().getPointsInPolygon(polygon, points);
        for(int i = 0; i < points.size(2); i++)
        {
            mInnerSamples.insert(points[i]);
        }

        mShell.getOuterTree().getPointsInPolygon(polygon, points);
        for(int i = 0; i < points.size(2); i++)
        {
            mOuterSamples.insert(points[i]);
        }
    }*/

    bool KernelRegion::contains(const vec3_t &point) const
    {
        vec4_t homo;
        homo[0] = point[0];
        homo[1] = point[1];
        homo[2] = point[2];
        homo[3] = 1;
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

    bool KernelRegion::isInvalidRegion(const vec3_t &point) const
    {
        for(int i = 0; i < mFaces.size(2); i++)
        {
            matrixr_t cell(3, 4);
            cell(colon(), 0) = mPoints(colon(), mFaces(0, i));
            cell(colon(), 1) = mPoints(colon(), mFaces(1, i));
            cell(colon(), 2) = mPoints(colon(), mFaces(2, i));
            cell(colon(), 3) = point;

            BaryComputer baryComputer(cell);
            std::set<size_t>::const_iterator itr = mInnerSamples.begin();
            for(int j = 0; j < mInnerSamples.size(); j++, ++itr)
            {
                vec4_t bary;
                baryComputer(mShell.mInnerShell(colon(), *itr), bary);
                if(min(bary) >= 0)
                {
                    double f0 = mTriangulation.getFValue(mFaces(0, i));
                    double f1 = mTriangulation.getFValue(mFaces(1, i));
                    double f2 = mTriangulation.getFValue(mFaces(2, i));
                    double f3 = mTriangulation.getFValue(mPointType);

                    double f = f0 * bary[0] + f1 * bary[1]+ f2 * bary[2] + f3 * bary[3];
                    if(f > 0)
                    {
                        return true;
                    }
                }
            }
            itr = mOuterSamples.begin();
            for(int j = 0; j < mOuterSamples.size(); j++, ++itr)
            {
                vec4_t bary;
                baryComputer(mShell.mOuterShell(colon(), *itr), bary);
                if(min(bary) >= 0)
                {
                    double f0 = mTriangulation.getFValue(mFaces(0, i));
                    double f1 = mTriangulation.getFValue(mFaces(1, i));
                    double f2 = mTriangulation.getFValue(mFaces(2, i));
                    double f3 = mTriangulation.getFValue(mPointType);

                    double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2] + f3 * bary[3];
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
