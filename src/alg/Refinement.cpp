#include "Refinement.h"
#include "BaryComputer.h"
#include <limits>
#include <wkylib/Geometry/geometry.h>
#include <zjucad/matrix/io.h>

using namespace zjucad::matrix;
using namespace WKYLIB::Geometry;

namespace SBV
{
    Refinement::Refinement(const Shell& shell, TriangulatedShell &output,
                           double alpha, double sampleRadius)
        : mShell(shell),
          mOutput(output),
          mAlpha(alpha),
          mSampleRadius(sampleRadius)
    {
        init();
    }

    void Refinement::init()
    {
        mOutput.vertType.clear();
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

        for(int i = 0; i < mShell.mOuterShell.size(2); i++)
        {
            if(xMax < mShell.mOuterShell(0, i))
            {
                xMax = mShell.mOuterShell(0, i);
            }
            if(xMin > mShell.mOuterShell(0, i))
            {
                xMin = mShell.mOuterShell(0, i);
            }
            if(yMax < mShell.mOuterShell(1, i))
            {
                yMax = mShell.mOuterShell(1, i);
            }
            if(yMin > mShell.mOuterShell(1, i))
            {
                yMin = mShell.mOuterShell(1, i);
            }
        }

        xMin *= 1.1;
        xMax *= 1.1;
        yMin *= 1.1;
        yMax *= 1.1;

        PointInfo info;
        info.pointType = PointType::POINT_BOUNDING_BOX;

        mDelaunay.insert(Point(xMax, yMax))->info() = info;
        mDelaunay.insert(Point(xMax, yMin))->info() = info;
        mDelaunay.insert(Point(xMin, yMax))->info() = info;
        mDelaunay.insert(Point(xMin, yMin))->info() = info;
    }

    void Refinement::initErrors()
    {
        mInnerError.resize(mShell.mInnerShell.size(2), 0);
        mOuterError.resize(mShell.mOuterShell.size(2), 0);
    }

    void Refinement::updateErrors()
    {
        double maxError = -std::numeric_limits<double>::max();
        for(auto iter = mDelaunay.finite_faces_begin(); iter != mDelaunay.finite_faces_end(); ++iter)
        {
            updatePointInCell(*iter);
            double maxErrorInCell = getError(iter->info().maxErrorPoint);
            if(maxError < maxErrorInCell)
            {
                maxError = maxErrorInCell;
                mNextInsertPoint = iter->info().maxErrorPoint;
            }
        }
    }

    void Refinement::updatePointInCell(Cell& cell)
    {
        FaceInfo& info = cell.info();
        if(isNewCell(cell) == false)
        {
            return;
        }

        info.v0 = cell.vertex(0)->info();
        info.v1 = cell.vertex(1)->info();
        info.v2 = cell.vertex(2)->info();

        double xmax, xmin, ymax, ymin;
        computeAABB(cell, xmax, xmin, ymax, ymin);
        matrixs_t sampleInner, sampleOuter;
        mShell.getInnerTree().getPointsInRange(xmin, xmax, ymin, ymax, sampleInner);
        mShell.getOuterTree().getPointsInRange(xmin, xmax, ymin, ymax, sampleOuter);

        double f0 = mOutput.getFValue(cell.vertex(0)->info().pointType);
        double f1 = mOutput.getFValue(cell.vertex(1)->info().pointType);
        double f2 = mOutput.getFValue(cell.vertex(2)->info().pointType);

        double maxError = std::numeric_limits<double>::lowest();
        PointInfo maxErrorPoint;

        const Point& p0 = cell.vertex(0)->point();
        const Point& p1 = cell.vertex(1)->point();
        const Point& p2 = cell.vertex(2)->point();
        matrixr_t triangle(2, 3);
        triangle(0, 0) = p0[0];
        triangle(1, 0) = p0[1];
        triangle(0, 1) = p1[0];
        triangle(1, 1) = p1[1];
        triangle(0, 2) = p2[0];
        triangle(1, 2) = p2[1];

        BaryComputer baryComputer(mShell, triangle, sampleInner, sampleOuter);
        matrixr_t barysInner, barysOuter;
        baryComputer.computeBary(barysInner, barysOuter);
        for(int i = 0; i < barysInner.size(2); i++)
        {
            const size_t index = sampleInner[i];
            const vec3_t bary = barysInner(colon(), i);

            if(min(bary) <= -1e-6)
                continue;

            double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2];
            double error = fabs(f + 1);
            if(maxError < error)
            {
                maxError = error;
                maxErrorPoint = PointInfo(POINT_INNER, index);
            }
            mInnerError[index] = error;
            info.points.push_back(PointInfo(POINT_INNER, index));
        }
        for(int i = 0; i < barysOuter.size(2); i++)
        {
            const size_t index = sampleOuter[i];
            const vec3_t bary = barysOuter(colon(), i);

            if(min(bary) <= -1e-6)
                continue;

            double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2];
            double error = fabs(f - 1);
            if(maxError < error)
            {
                maxError = error;
                maxErrorPoint = PointInfo(POINT_OUTER, index);
            }
            mOuterError[index] = error;
            info.points.push_back(PointInfo(POINT_OUTER, index));
        }

        info.maxErrorPoint = maxErrorPoint;
    }

    void Refinement::computeAABB(const Cell &cell, double &xmax, double &xmin, double &ymax, double &ymin)
    {
        const Point& p0 = cell.vertex(0)->point();
        const Point& p1 = cell.vertex(1)->point();
        const Point& p2 = cell.vertex(2)->point();

        xmax = std::numeric_limits<double>::lowest();
        xmin = std::numeric_limits<double>::max();
        ymax = std::numeric_limits<double>::lowest();
        ymin = std::numeric_limits<double>::max();

        xmax = std::max(xmax, std::max(std::max(p0[0], p1[0]), p2[0]));
        xmin = std::min(xmin, std::min(std::min(p0[0], p1[0]), p2[0]));
        ymax = std::max(ymax, std::max(std::max(p0[1], p1[1]), p2[1]));
        ymin = std::min(ymin, std::min(std::min(p0[1], p1[1]), p2[1]));
    }

    double Refinement::computeFValue(const vec2_t &point, const Cell &cell)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);

        Point& p0 = cell.vertex(0)->point();
        Point& p1 = cell.vertex(1)->point();
        Point& p2 = cell.vertex(2)->point();

        vec2_t a, b, c;
        a[0] = p0[0];
        a[1] = p0[1];
        b[0] = p1[0];
        b[1] = p1[1];
        c[0] = p2[0];
        c[1] = p2[1];

        vec3_t bary;
        if(barycentric_2D(a, b, c, point, bary))
        {
            //the point is inside the tetrahedron
            double f0 = getFValue(vh0);
            double f1 = getFValue(vh1);
            double f2 = getFValue(vh2);

            double ans = f0 * bary[0] + f1 * bary[1] + f2 * bary[2];
            //assert(ans < 1 || fabs(ans-1) < 1e-3);
            //assert(ans > -1 || fabs(ans+1) < 1e-3);
            return ans;
        }

        //the point is outside the tetrahedron, so returns an error
        return std::numeric_limits<double>::max();
    }

    double Refinement::getFValue(const VertexHandle& vh)
    {
        const PointInfo& info = vh->info();
        return mOutput.getFValue(info.pointType);
    }

    double Refinement::getError(const PointInfo &point)
    {
        switch(point.pointType)
        {
        case PointType::POINT_INNER:
            return mInnerError[point.index];
        case PointType::POINT_OUTER:
            return mOuterError[point.index];
        case PointType::POINT_UNKNOWN:
            return std::numeric_limits<double>::max();
        default:
            //this should not be run, otherwise there exists logic error
            throw std::logic_error("In getError(), pointType should be POINT_INNER or POINT_OUTER");
        }
    }

    void Refinement::getPointMatrix(const PointInfo &point, vec2_t &pointMatrix)
    {
        switch(point.pointType)
        {
        case PointType::POINT_INNER:
            pointMatrix = mShell.mInnerShell(colon(), point.index);
            break;
        case PointType::POINT_OUTER:
            pointMatrix = mShell.mOuterShell(colon(), point.index);
            break;
        default:
            //this should not be runned, otherwise there exists logic error
            throw std::logic_error("In getPointMatrix(), pointType should be POINT_INNER or POINT_OUTER");
        }
    }

    bool Refinement::refine()
    {
        int iterCount = 0;
        while(!isFinished() && iterCount <= (mShell.mInnerShell.size(2) + mShell.mOuterShell.size(2)))
        {
            iterCount++;
            std::cout << "Iteration " << iterCount << " ..." << std::endl;

            vec2_t point;
            getPointMatrix(mNextInsertPoint, point);
            mDelaunay.insert(Point(point[0], point[1]))->info() = mNextInsertPoint;
            updateErrors();
        }

        std::cout << "Finished refinement after " << iterCount << " iterations." << std::endl;

        //organize output data
        mOutput.vertices.resize(2, mDelaunay.number_of_vertices());
        mOutput.triangles.resize(3, mDelaunay.number_of_faces());
        mOutput.vertType.reserve(mDelaunay.number_of_vertices());
        int i = 0;
        for(auto iter = mDelaunay.finite_vertices_begin(); iter != mDelaunay.finite_vertices_end(); ++iter, ++i)
        {
            mOutput.vertices(0, i) = iter->point()[0];
            mOutput.vertices(1, i) = iter->point()[1];

            PointInfo& info = iter->info();
            iter->info().indexInDelaunay = i;
            mOutput.vertType.push_back(info.pointType);
        }
        i = 0;
        for(auto iter = mDelaunay.finite_faces_begin(); iter != mDelaunay.finite_faces_end(); ++iter, ++i)
        {
            const PointInfo& info0 = iter->vertex(0)->info();
            const PointInfo& info1 = iter->vertex(1)->info();
            const PointInfo& info2 = iter->vertex(2)->info();

            mOutput.triangles(0, i) = info0.indexInDelaunay;
            mOutput.triangles(1, i) = info1.indexInDelaunay;
            mOutput.triangles(2, i) = info2.indexInDelaunay;
        }

        return true;
    }

    bool Refinement::isFinished()
    {
        //check for condition 1.
        for(int i = 0; i < mInnerError.size(); i++)
        {
            if(mInnerError[i] > 1 - mAlpha)
            {
                return false;
            }
        }

        for(int i = 0; i < mOuterError.size(); i++)
        {
            if(mOuterError[i] > 1 - mAlpha)
            {
                return false;
            }
        }

        //check for condition 2 and 3.
        for(auto iter = mDelaunay.finite_faces_begin(); iter != mDelaunay.finite_faces_end(); ++iter)
        {
            const Cell& cell = *iter;
            const VertexHandle& vh0 = cell.vertex(0);
            const VertexHandle& vh1 = cell.vertex(1);
            const VertexHandle& vh2 = cell.vertex(2);

            if(getFValue(vh0) == getFValue(vh1) && getFValue(vh0) == getFValue(vh2))
            {
                continue;
            }

            /*
            //condition 2
            if(computeHeight(cell) < 2 * mSampleRadius / mAlpha)
            {
                return false;
            }*/

            if(checkCondition3(cell) == false)
            {
                return false;
            }
        }

        return true;
    }

    bool Refinement::isNewCell(const Cell &cell)
    {
        const FaceInfo& info = cell.info();
        const PointInfo& p0 = cell.vertex(0)->info();
        const PointInfo& p1 = cell.vertex(1)->info();
        const PointInfo& p2 = cell.vertex(2)->info();

        return !((p0 == info.v0 || p0 == info.v1 || p0 == info.v2)
                &&(p1 == info.v0 || p1 == info.v1 || p1 == info.v2)
                &&(p2 == info.v0 || p2 == info.v1 || p2 == info.v2));
    }

    double Refinement::computeHeight(const Cell &cell)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);
        const Point& p0 = vh0->point();
        const Point& p1 = vh1->point();
        const Point& p2 = vh2->point();

        vec3_t a, b, c;

        if(getFValue(vh0) == getFValue(vh1))
        {
            a[0] = p2[0];
            a[1] = p2[1];
            a[2] = 1;
            b[0] = p0[0];
            b[1] = p0[1];
            b[2] = 1;
            c[0] = p1[0];
            c[1] = p1[1];
            c[2] = 1;
        }
        else if(getFValue(vh0) == getFValue(vh2))
        {
            a[0] = p1[0];
            a[1] = p1[1];
            a[2] = 1;
            b[0] = p0[0];
            b[1] = p0[1];
            b[2] = 1;
            c[0] = p2[0];
            c[1] = p2[1];
            c[2] = 1;
        }
        else if(getFValue(vh1) == getFValue(vh2))
        {
            a[0] = p0[0];
            a[1] = p0[1];
            a[2] = 1;
            b[0] = p1[0];
            b[1] = p1[1];
            b[2] = 1;
            c[0] = p2[0];
            c[1] = p2[1];
            c[2] = 1;
        }

        //compute height from a to bc.
        vec3_t ba = a - b;
        vec3_t bc = c - b;
        return norm(cross(ba, bc)) / norm(bc);
    }

    bool Refinement::checkCondition3(const Cell &cell)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);
        const Point& p0 = vh0->point();
        const Point& p1 = vh1->point();
        const Point& p2 = vh2->point();

        vec2_t a, b, c;
        a[0] = p0[0];
        a[1] = p0[1];
        b[0] = p1[0];
        b[1] = p1[1];
        c[0] = p2[0];
        c[1] = p2[1];

        vec2_t center = (a + b + c) / 3.0;
        constexpr double k = sqrt(0.7);
        vec2_t newA = k * a + (1 - k) * center;
        vec2_t newB = k * b + (1 - k) * center;
        vec2_t newC = k * c + (1 - k) * center;

        const KdTreeWrap& innerTree = mShell.getInnerTree();
        size_t nearA, nearB, nearC;
        nearA = innerTree.getNearestPoint(newA);
        nearB = innerTree.getNearestPoint(newB);
        nearC = innerTree.getNearestPoint(newC);
        if(checkClassification(cell, mShell.mInnerShell(colon(), nearA), false) == false)
        {
            return false;
        }
        if(checkClassification(cell, mShell.mInnerShell(colon(), nearB), false) == false)
        {
            return false;
        }
        if(checkClassification(cell, mShell.mInnerShell(colon(), nearC), false) == false)
        {
            return false;
        }

        const KdTreeWrap& outerTree = mShell.getOuterTree();
        nearA = outerTree.getNearestPoint(newA);
        nearB = outerTree.getNearestPoint(newB);
        nearC = outerTree.getNearestPoint(newC);
        if(checkClassification(cell, mShell.mOuterShell(colon(), nearA), true) == false)
        {
            return false;
        }
        if(checkClassification(cell, mShell.mOuterShell(colon(), nearB), true) == false)
        {
            return false;
        }
        if(checkClassification(cell, mShell.mOuterShell(colon(), nearC), true) == false)
        {
            return false;
        }

        return true;
    }

    bool Refinement::checkClassification(const Cell &cell, const vec2_t &point, bool isOuter)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);
        const Point& p0 = vh0->point();
        const Point& p1 = vh1->point();
        const Point& p2 = vh2->point();

        vec2_t a, b, c;
        a[0] = p0[0];
        a[1] = p0[1];
        b[0] = p1[0];
        b[1] = p1[1];
        c[0] = p2[0];
        c[1] = p2[1];

        vec3_t bary;
        barycentric_2D(a, b, c, point, bary);
        double fvalue = getFValue(vh0) * bary[0] + getFValue(vh1) * bary[1] + getFValue(vh2) * bary[2];
        bool result;
        if(isOuter)
        {
            result = fvalue > mAlpha;
        }
        else
        {
            result = fvalue < -mAlpha;
        }
        return result;
    }
}
