#include "Refinement.h"
#include <limits>
#include <wkylib/geometry.h>

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
        double xMax = std::numeric_limits<double>::lowest();
        double xMin = std::numeric_limits<double>::max();
        double yMax = std::numeric_limits<double>::lowest();
        double yMin = std::numeric_limits<double>::max();
#ifndef VER_2D
        double zMax = std::numeric_limits<double>::lowest();
        double zMin = std::numeric_limits<double>::max();
#endif

        for(size_t i = 0; i < mShell.mOuterShell.size(); i++)
        {
            if(xMax < mShell.mOuterShell[i][0])
            {
                xMax = mShell.mOuterShell[i][0];
            }
            if(xMin > mShell.mOuterShell[i][0])
            {
                xMin = mShell.mOuterShell[i][0];
            }
            if(yMax < mShell.mOuterShell[i][1])
            {
                yMax = mShell.mOuterShell[i][1];
            }
            if(yMin > mShell.mOuterShell[i][1])
            {
                yMin = mShell.mOuterShell[i][1];
            }
#ifndef VER_2D
            if(zMax < mShell.mOuterShell[i][2])
            {
                zMax = mShell.mOuterShell[i][2];
            }
            if(zMin > mShell.mOuterShell[i][2])
            {
                zMin = mShell.mOuterShell[i][2];
            }
#endif
        }

        if(xMin != std::numeric_limits<double>::max())
            xMin *= 1.1;
        if(xMax != std::numeric_limits<double>::lowest())
            xMax *= 1.1;
        if(yMin != std::numeric_limits<double>::max())
            yMin *= 1.1;
        if(yMax != std::numeric_limits<double>::lowest())
            yMax *= 1.1;
#ifndef VER_2D
        if(zMin != std::numeric_limits<double>::max())
            zMin *= 1.1;
        if(zMax != std::numeric_limits<double>::lowest())
            zMax *= 1.1;
#endif

        PointInfo info;
        info.pointType = PointType::POINT_BOUNDING_BOX;

#ifdef VER_2D
        mDelaunay.insert(DPoint(xMax, yMax))->info() = info;
        mDelaunay.insert(DPoint(xMax, yMin))->info() = info;
        mDelaunay.insert(DPoint(xMin, yMax))->info() = info;
        mDelaunay.insert(DPoint(xMin, yMin))->info() = info;
#else
        mDelaunay.insert(DPoint(xMax, yMax, zMax))->info() = info;
        mDelaunay.insert(DPoint(xMax, yMax, zMin))->info() = info;
        mDelaunay.insert(DPoint(xMax, yMin, zMax))->info() = info;
        mDelaunay.insert(DPoint(xMax, yMin, zMin))->info() = info;
        mDelaunay.insert(DPoint(xMin, yMax, zMax))->info() = info;
        mDelaunay.insert(DPoint(xMin, yMax, zMin))->info() = info;
        mDelaunay.insert(DPoint(xMin, yMin, zMax))->info() = info;
        mDelaunay.insert(DPoint(xMax, yMin, zMin))->info() = info;
#endif
    }

    void Refinement::initErrors()
    {
        mInnerError.resize(mShell.mInnerShell.size(), 0);
        mOuterError.resize(mShell.mOuterShell.size(), 0);
    }

    void Refinement::updateErrors()
    {
        double maxError = std::numeric_limits<double>::lowest();
#ifdef VER_2D
        for(auto iter = mDelaunay.finite_faces_begin(); iter != mDelaunay.finite_faces_end(); ++iter)
#else
        for(auto iter = mDelaunay.finite_cells_begin(); iter != mDelaunay.finite_cells_end(); ++iter)
#endif
        {
            Cell& cell = *iter;
            updatePointInCell(cell);
            double maxErrorInCell = getError(cell.info().maxErrorPoint);
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

#ifdef VER_2D
        double xmax, xmin, ymax, ymin;
        computeAABB(cell, xmax, xmin, ymax, ymin);
        std::vector<size_t> sampleInner, sampleOuter;
        mShell.getInnerTree().getPointsInRange(xmin, xmax, ymin, ymax, sampleInner);
        mShell.getOuterTree().getPointsInRange(xmin, xmax, ymin, ymax, sampleOuter);
#else
        double xmax, xmin, ymax, ymin, zmax, zmin;
        computeAABB(cell, xmax, xmin, ymax, ymin, zmax, zmin);
        std::vector<size_t> sampleInner, sampleOuter;
        mShell.getInnerTree().getPointsInRange(xmin, xmax, ymin, ymax, zmin, zmax, sampleInner);
        mShell.getOuterTree().getPointsInRange(xmin, xmax, ymin, ymax, zmin, zmax, sampleOuter);
#endif

        double maxError = std::numeric_limits<double>::lowest();
        PointInfo maxErrorPoint;
        for(size_t i = 0; i < sampleInner.size(); i++)
        {
            const size_t index = sampleInner[i];
            double f = computeFValue(mShell.mInnerShell[index], cell);
            if(f == std::numeric_limits<double>::max())
            {
                //the point is outside the tetrahedron
                continue;
            }
            double error = fabs(f + 1);
            if(maxError < error)
            {
                maxError = error;
                maxErrorPoint = PointInfo(POINT_INNER, index);
            }
            mInnerError[index] = error;
            info.points.push_back(PointInfo(POINT_INNER, index));
        }

        for(size_t i = 0; i < sampleOuter.size(); i++)
        {
            const size_t index = sampleOuter[i];
            double f = computeFValue(mShell.mOuterShell[index], cell);
            if(f == std::numeric_limits<double>::max())
            {
                //the point is outside the tetrahedron
                continue;
            }
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

#ifdef VER_2D
    void Refinement::computeAABB(const Cell &cell, double &xmax, double &xmin, double &ymax, double &ymin)
    {
        const DPoint& p0 = cell.vertex(0)->point();
        const DPoint& p1 = cell.vertex(1)->point();
        const DPoint& p2 = cell.vertex(2)->point();

        xmax = std::numeric_limits<double>::lowest();
        xmin = std::numeric_limits<double>::max();
        ymax = std::numeric_limits<double>::lowest();
        ymin = std::numeric_limits<double>::max();

        xmax = std::max(xmax, std::max(std::max(p0[0], p1[0]), p2[0]));
        xmin = std::min(xmin, std::min(std::min(p0[0], p1[0]), p2[0]));
        ymax = std::max(ymax, std::max(std::max(p0[1], p1[1]), p2[1]));
        ymin = std::min(ymin, std::min(std::min(p0[1], p1[1]), p2[1]));
    }   
#else
    void Refinement::computeAABB(const Cell &cell, double &xmax, double &xmin, double &ymax, double &ymin, double& zmax, double& zmin)
    {
        const DPoint& p0 = cell.vertex(0)->point();
        const DPoint& p1 = cell.vertex(1)->point();
        const DPoint& p2 = cell.vertex(2)->point();
        const DPoint& p3 = cell.vertex(3)->point();

        xmax = std::numeric_limits<double>::lowest();
        xmin = std::numeric_limits<double>::max();
        ymax = std::numeric_limits<double>::lowest();
        ymin = std::numeric_limits<double>::max();
        zmax = std::numeric_limits<double>::lowest();
        zmin = std::numeric_limits<double>::max();

        xmax = std::max(xmax, std::max(std::max(std::max(p0[0], p1[0]), p2[0]), p3[0]));
        xmin = std::min(xmin, std::min(std::min(std::min(p0[0], p1[0]), p2[0]), p3[0]));
        ymax = std::max(ymax, std::max(std::max(std::max(p0[1], p1[1]), p2[1]), p3[1]));
        ymin = std::min(ymin, std::min(std::min(std::min(p0[1], p1[1]), p2[1]), p3[1]));
        zmax = std::max(zmax, std::max(std::max(std::max(p0[2], p1[2]), p2[2]), p3[2]));
        zmin = std::min(zmin, std::min(std::min(std::min(p0[2], p1[2]), p2[2]), p3[2]));
    }
#endif

    //optimize this function causes crash
#ifdef VER_2D
    double Refinement::computeFValue(const Point &point, const Cell &cell)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);

        DPoint& p0 = cell.vertex(0)->point();
        DPoint& p1 = cell.vertex(1)->point();
        DPoint& p2 = cell.vertex(2)->point();

        Point a, b, c;
        a[0] = p0[0];
        a[1] = p0[1];
        b[0] = p1[0];
        b[1] = p1[1];
        c[0] = p2[0];
        c[1] = p2[1];

        Eigen::Vector3d bary;
        if(WKYLIB::barycentric_2D(a, b, c, point, bary))
        {
            //the point is inside the tetrahedron
            double f0 = getFValue(vh0);
            double f1 = getFValue(vh1);
            double f2 = getFValue(vh2);

            double ans = f0 * bary[0] + f1 * bary[1] + f2 * bary[2];
            return ans;
        }

        //the point is outside the tetrahedron, so returns an error
        return std::numeric_limits<double>::max();
    }
#else
    double Refinement::computeFValue(const Point &point, const Cell &cell)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);
        const VertexHandle& vh3 = cell.vertex(3);

        DPoint& p0 = cell.vertex(0)->point();
        DPoint& p1 = cell.vertex(1)->point();
        DPoint& p2 = cell.vertex(2)->point();
        DPoint& p3 = cell.vertex(3)->point();

        Point a, b, c, d;
        a[0] = p0[0];
        a[1] = p0[1];
        a[2] = p0[2];
        b[0] = p1[0];
        b[1] = p1[1];
        b[2] = p1[2];
        c[0] = p2[0];
        c[1] = p2[1];
        c[2] = p2[2];
        d[0] = p3[0];
        d[1] = p3[1];
        d[2] = p3[2];

        Eigen::Vector4d bary;
        if(WKYLIB::barycentric_tetra(a, b, c, d, point, bary))
        {
            //the point is inside the tetrahedron
            double f0 = getFValue(vh0);
            double f1 = getFValue(vh1);
            double f2 = getFValue(vh2);
            double f3 = getFValue(vh3);

            double ans = f0 * bary[0] + f1 * bary[1] + f2 * bary[2] + f3 * bary[3];
            return ans;
        }

        //the point is outside the tetrahedron, so returns an error
        return std::numeric_limits<double>::max();
    }
#endif

    //optimize this function causes crash
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
        default:
            //this should not be run, otherwise there exists logic error
            throw std::logic_error("In getError(), pointType should be POINT_INNER or POINT_OUTER");
        }
    }

    void Refinement::getPointMatrix(const PointInfo &point, Point &pointMatrix)
    {
        switch(point.pointType)
        {
        case PointType::POINT_INNER:
            pointMatrix = mShell.mInnerShell[point.index];
            break;
        case PointType::POINT_OUTER:
            pointMatrix = mShell.mOuterShell[point.index];
            break;
        default:
            //this should not be runned, otherwise there exists logic error
            throw std::logic_error("In getPointMatrix(), pointType should be POINT_INNER or POINT_OUTER");
        }
    }

    bool Refinement::refine()
    {
        size_t iterCount = 0;
        while(!isFinished() && iterCount <= (mShell.mInnerShell.size() + mShell.mOuterShell.size()))
        {
            iterCount++;
            std::cout << "Iteration " << iterCount << " ..." << std::endl;

            Point point;
            getPointMatrix(mNextInsertPoint, point);
#ifdef VER_2D
            mDelaunay.insert(DPoint(point[0], point[1]))->info() = mNextInsertPoint;
#else
            mDelaunay.insert(DPoint(point[0], point[1], point[2]))->info() = mNextInsertPoint;
#endif
            updateErrors();
        }

        std::cout << "Finished refinement after " << iterCount << " iterations." << std::endl;

        //organize output data
        mOutput.vertices.resize(mDelaunay.number_of_vertices());
#ifdef VER_2D
        mOutput.triangles.resize(mDelaunay.number_of_faces());
#else
        mOutput.triangles.resize(mDelaunay.number_of_cells());
#endif
        mOutput.vertType.reserve(mDelaunay.number_of_vertices());
        int i = 0;
        for(auto iter = mDelaunay.finite_vertices_begin(); iter != mDelaunay.finite_vertices_end(); ++iter, ++i)
        {
            mOutput.vertices[i][0] = iter->point()[0];
            mOutput.vertices[i][1] = iter->point()[1];
#ifndef VER_2D
            mOutput.vertices[i][2] = iter->point()[2];
#endif

            PointInfo& info = iter->info();
            iter->info().indexInDelaunay = i;
            mOutput.vertType.push_back(info.pointType);
        }
        i = 0;
#ifdef VER_2D
        for(auto iter = mDelaunay.finite_faces_begin(); iter != mDelaunay.finite_faces_end(); ++iter, ++i)
        {
            const PointInfo& info0 = iter->vertex(0)->info();
            const PointInfo& info1 = iter->vertex(1)->info();
            const PointInfo& info2 = iter->vertex(2)->info();

            mOutput.triangles[i][0] = info0.indexInDelaunay;
            mOutput.triangles[i][1] = info1.indexInDelaunay;
            mOutput.triangles[i][2] = info2.indexInDelaunay;
        }
#else
        for(auto iter = mDelaunay.finite_cells_begin(); iter != mDelaunay.finite_cells_end(); ++iter, ++i)
        {
            const PointInfo& info0 = iter->vertex(0)->info();
            const PointInfo& info1 = iter->vertex(1)->info();
            const PointInfo& info2 = iter->vertex(2)->info();
            const PointInfo& info3 = iter->vertex(3)->info();

            mOutput.triangles[i][0] = info0.indexInDelaunay;
            mOutput.triangles[i][1] = info1.indexInDelaunay;
            mOutput.triangles[i][2] = info2.indexInDelaunay;
            mOutput.triangles[i][3] = info3.indexInDelaunay;
        }
#endif

        return true;
    }

    bool Refinement::isFinished()
    {
        //check for condition 1.
        for(size_t i = 0; i < mInnerError.size(); i++)
        {
            if(mInnerError[i] > 1 - mAlpha)
            {
                return false;
            }
        }

        for(size_t i = 0; i < mOuterError.size(); i++)
        {
            if(mOuterError[i] > 1 - mAlpha)
            {
                return false;
            }
        }

        //check for condition 2 and 3.
#ifdef VER_2D
        for(auto iter = mDelaunay.finite_faces_begin(); iter != mDelaunay.finite_faces_end(); ++iter)
#else
        for(auto iter = mDelaunay.finite_cells_begin(); iter != mDelaunay.finite_cells_end(); ++iter)
#endif
        {
            const Cell& cell = *iter;
            const VertexHandle& vh0 = cell.vertex(0);
            const VertexHandle& vh1 = cell.vertex(1);
            const VertexHandle& vh2 = cell.vertex(2);
#ifndef VER_2D
            const VertexHandle& vh3 = cell.vertex(2);
#endif

#ifdef VER_2D
            if(getFValue(vh0) == getFValue(vh1) && getFValue(vh0) == getFValue(vh2))
            {
                continue;
            }
#else
            if(getFValue(vh0) == getFValue(vh1) && getFValue(vh0) == getFValue(vh2) && getFValue(vh0) == getFValue(vh3))
            {
                continue;
            }
#endif

            //condition 2
            //if(computeHeight(cell) < 2 * mSampleRadius / mAlpha)
            //{
            //    return false;
            //}

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
        const DPoint& p0 = vh0->point();
        const DPoint& p1 = vh1->point();
        const DPoint& p2 = vh2->point();

        Eigen::Vector3d a, b, c;

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
        Eigen::Vector3d ba = a - b;
        Eigen::Vector3d bc = c - b;
        return (ba.cross(bc)).norm() / bc.norm();
    }

    bool Refinement::checkCondition3(const Cell &cell)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);
        const DPoint& p0 = vh0->point();
        const DPoint& p1 = vh1->point();
        const DPoint& p2 = vh2->point();

        Point a, b, c;
        a[0] = p0[0];
        a[1] = p0[1];
        b[0] = p1[0];
        b[1] = p1[1];
        c[0] = p2[0];
        c[1] = p2[1];

        Point center = (a + b + c) / 3.0;
        constexpr double k = sqrt(0.7);
        Point newA = k * a + (1 - k) * center;
        Point newB = k * b + (1 - k) * center;
        Point newC = k * c + (1 - k) * center;

        const KdTreeWrap& innerTree = mShell.getInnerTree();
        size_t nearA, nearB, nearC;
        nearA = innerTree.getNearestPoint(newA);
        nearB = innerTree.getNearestPoint(newB);
        nearC = innerTree.getNearestPoint(newC);
        if(checkClassification(cell, mShell.mInnerShell[nearA], false) == false)
        {
            return false;
        }
        if(checkClassification(cell, mShell.mInnerShell[nearB], false) == false)
        {
            return false;
        }
        if(checkClassification(cell, mShell.mInnerShell[nearC], false) == false)
        {
            return false;
        }

        const KdTreeWrap& outerTree = mShell.getOuterTree();
        nearA = outerTree.getNearestPoint(newA);
        nearB = outerTree.getNearestPoint(newB);
        nearC = outerTree.getNearestPoint(newC);
        if(checkClassification(cell, mShell.mOuterShell[nearA], true) == false)
        {
            return false;
        }
        if(checkClassification(cell, mShell.mOuterShell[nearB], true) == false)
        {
            return false;
        }
        if(checkClassification(cell, mShell.mOuterShell[nearC], true) == false)
        {
            return false;
        }

        return true;
    }

#ifdef VER_2D
    bool Refinement::checkClassification(const Cell &cell, const Point &point, bool isOuter)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);
        const DPoint& p0 = vh0->point();
        const DPoint& p1 = vh1->point();
        const DPoint& p2 = vh2->point();

        Point a, b, c;
        a[0] = p0[0];
        a[1] = p0[1];
        b[0] = p1[0];
        b[1] = p1[1];
        c[0] = p2[0];
        c[1] = p2[1];

        Eigen::Vector3d bary;
        WKYLIB::barycentric_2D(a, b, c, point, bary);
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
#else
    bool Refinement::checkClassification(const Cell &cell, const Point &point, bool isOuter)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);
        const VertexHandle& vh3 = cell.vertex(3);

        DPoint& p0 = cell.vertex(0)->point();
        DPoint& p1 = cell.vertex(1)->point();
        DPoint& p2 = cell.vertex(2)->point();
        DPoint& p3 = cell.vertex(3)->point();

        Point a, b, c, d;
        a[0] = p0[0];
        a[1] = p0[1];
        a[2] = p0[2];
        b[0] = p1[0];
        b[1] = p1[1];
        b[2] = p1[2];
        c[0] = p2[0];
        c[1] = p2[1];
        c[2] = p2[2];
        d[0] = p3[0];
        d[1] = p3[1];
        d[2] = p3[2];

        Eigen::Vector4d bary;
        WKYLIB::barycentric_tetra(a, b, c, d, point, bary);
        double fvalue = getFValue(vh0) * bary[0] + getFValue(vh1) * bary[1] + getFValue(vh2) * bary[2] + getFValue(vh3) * bary[3];
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
#endif
}



