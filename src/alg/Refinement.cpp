#include "Refinement.h"
#include <limits>
#include <wkylib/geometry.h>

namespace SBV
{
    Refinement::Refinement(const matrixr_t &innerShell, const matrixr_t &outerShell)
        : mInnerShell(innerShell),
          mOuterShell(outerShell)
    {
        init();
    }

    void Refinement::init()
    {
        computeBoundingBox();
        initErrors();
        updateErrors();
    }

    void Refinement::computeBoundingBox()
    {
        double xMax = std::numeric_limits<double>::min();
        double xMin = std::numeric_limits<double>::max();
        double yMax = std::numeric_limits<double>::min();
        double yMin = std::numeric_limits<double>::max();
        double zMax = std::numeric_limits<double>::min();
        double zMin = std::numeric_limits<double>::max();

        for(int i = 0; i < mOuterShell.size(2); i++)
        {
            if(xMax < mOuterShell(0, i))
            {
                xMax = mOuterShell(0, i);
            }
            if(xMin > mOuterShell(0, i))
            {
                xMin = mOuterShell(0, i);
            }
            if(yMax < mOuterShell(1, i))
            {
                yMax = mOuterShell(1, i);
            }
            if(yMin > mOuterShell(1, i))
            {
                yMin = mOuterShell(1, i);
            }
            if(zMax < mOuterShell(2, i))
            {
                zMax = mOuterShell(2, i);
            }
            if(zMin > mOuterShell(2, i))
            {
                zMin = mOuterShell(2, i);
            }
        }

        PointInfo info;
        info.pointType = PointType::POINT_BOUNDING_BOX;

        mDelaunay.insert(Point(xMax, yMax, zMax))->info() = info;
        mDelaunay.insert(Point(xMax, yMax, zMin))->info() = info;
        mDelaunay.insert(Point(xMax, yMin, zMax))->info() = info;
        mDelaunay.insert(Point(xMax, yMin, zMin))->info() = info;
        mDelaunay.insert(Point(xMin, yMax, zMax))->info() = info;
        mDelaunay.insert(Point(xMin, yMax, zMin))->info() = info;
        mDelaunay.insert(Point(xMin, yMin, zMax))->info() = info;
        mDelaunay.insert(Point(xMin, yMin, zMin))->info() = info;
    }

    void Refinement::initErrors()
    {
        mInnerError.reserve(mInnerShell.size(2));
        for(int i = 0; i < mInnerShell.size(2); i++)
        {
            mInnerError.push_back(0);
        }

        mOuterError.reserve(mOuterShell.size(2));
        for(int i = 0; i < mOuterShell.size(2); i++)
        {
            mOuterError.push_back(0);
        }
    }

    void Refinement::updateErrors()
    {
        for(auto iter = mDelaunay.finite_cells_begin(); iter != mDelaunay.finite_cells_end(); ++iter)
        {
            updatePointInCell(*iter);
        }
    }

    void Refinement::updatePointInCell(const Cell& cell)
    {
        for(int i = 0; i < mInnerShell.size(2); i++)
        {
            double f = computeFValue(mInnerShell(colon(), i), cell);
            if(f == std::numeric_limits<double>::max())
            {
                //the point is outside the tetrahedron
                continue;
            }
            mInnerError[i] = fabs(f + 1);
        }

        for(int i = 0; i < mOuterShell.size(2); i++)
        {
            double f = computeFValue(mOuterShell(colon(), i), cell);
            if(f == std::numeric_limits<double>::max())
            {
                //the point is outside the tetrahedron
                continue;
            }
            mOuterError[i] = fabs(f - 1);
        }
    }

    double Refinement::computeFValue(const matrixr_t &point, const Cell &cell)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);
        const VertexHandle& vh3 = cell.vertex(3);

        Point& p0 = cell.vertex(0)->point();
        Point& p1 = cell.vertex(1)->point();
        Point& p2 = cell.vertex(2)->point();
        Point& p3 = cell.vertex(3)->point();

        matrixr_t tetra(3, 4);
        tetra(0, 0) = p0[0];
        tetra(1, 0) = p0[1];
        tetra(2, 0) = p0[2];
        tetra(0, 1) = p1[0];
        tetra(1, 1) = p1[1];
        tetra(2, 1) = p1[2];
        tetra(0, 2) = p2[0];
        tetra(1, 2) = p2[1];
        tetra(2, 2) = p2[2];
        tetra(0, 3) = p3[0];
        tetra(1, 3) = p3[1];
        tetra(2, 3) = p3[2];

        matrixr_t bary;
        if(WKYLIB::barycentric_tetra(point, tetra, bary))
        {
            //the point is inside the tetrahedron
            double f0 = getFValue(vh0);
            double f1 = getFValue(vh1);
            double f2 = getFValue(vh2);
            double f3 = getFValue(vh3);

            return f0 * bary[0] + f1 * bary[1] + f2 * bary[2] + f3 * bary[3];
        }

        //the point is outside the tetrahedron, so returns an error
        return std::numeric_limits<double>::max();
    }

    double Refinement::getFValue(const VertexHandle& vh)
    {
        const PointInfo& info = vh->info();
        switch(info.pointType)
        {
        case PointType::POINT_BOUNDING_BOX:
        case PointType::POINT_OUTER:
            return 1;
        case PointType::POINT_INNER:
            return -1;
        }
    }

    bool Refinement::refine(std::vector<size_t> &output_refinement)
    {

        while(!isFinished())
        {

        }

        return true;
    }

    bool Refinement::isFinished()
    {
        for(int i = 0; i < mInnerError.size(); i++)
        {
            if(mInnerError[i] > 0.8)
            {
                return false;
            }
        }

        for(int i = 0; i < mOuterError.size(); i++)
        {
            if(mOuterError[i] > 0.8)
            {
                return false;
            }
        }

        return true;
    }
}
