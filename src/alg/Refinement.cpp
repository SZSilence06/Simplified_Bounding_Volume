#include "Refinement.h"
#include "BaryComputer.h"
#include "NormalChecker.h"
#include <limits>
#include <wkylib/geometry.h>
#include <zjucad/matrix/io.h>
#include <CGAL/Kernel/global_functions_3.h>

int bad_cell = -1;
bool is_outer = false;
int bad_sample_index = -1;

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
        double zMax = std::numeric_limits<double>::lowest();
        double zMin = std::numeric_limits<double>::max();

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
            if(zMax < mShell.mOuterShell(2, i))
            {
                zMax = mShell.mOuterShell(2, i);
            }
            if(zMin > mShell.mOuterShell(2, i))
            {
                zMin = mShell.mOuterShell(2, i);
            }
        }

        xMin *= 1.1;
        xMax *= 1.1;
        yMin *= 1.1;
        yMax *= 1.1;
        zMin *= 1.1;
        zMax *= 1.1;

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
        mInnerError.resize(mShell.mInnerShell.size(2), 0);
        mOuterError.resize(mShell.mOuterShell.size(2), 0);
        mInnerExists.resize(mShell.mInnerShell.size(2), false);
        mOuterExists.resize(mShell.mOuterShell.size(2), false);
    }

    void Refinement::updateErrors()
    {
        double maxError = std::numeric_limits<double>::lowest();
        for(auto iter = mDelaunay.finite_cells_begin(); iter != mDelaunay.finite_cells_end(); ++iter)
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
        CellInfo& info = cell.info();
        if(info.isNew == false)
        {
            return;
        }

        info.isNew = false;

        double xmax, xmin, ymax, ymin, zmax, zmin;
        computeAABB(cell, xmax, xmin, ymax, ymin, zmax, zmin);
        matrixs_t sampleInner, sampleOuter;
        mShell.getInnerTree().getPointsInRange(xmin, xmax, ymin, ymax, zmin, zmax, sampleInner);
        mShell.getOuterTree().getPointsInRange(xmin, xmax, ymin, ymax, zmin, zmax, sampleOuter);

        double f0 = mOutput.getFValue(cell.vertex(0)->info().pointType);
        double f1 = mOutput.getFValue(cell.vertex(1)->info().pointType);
        double f2 = mOutput.getFValue(cell.vertex(2)->info().pointType);
        double f3 = mOutput.getFValue(cell.vertex(3)->info().pointType);

        double maxError = std::numeric_limits<double>::lowest();
        PointInfo maxErrorPoint;

        matrixr_t tetra(3, 4);
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                tetra(j, i) = cell.vertex(i)->point()[j];
            }
        }

        BaryComputer baryComputer(tetra);
        for(int i = 0; i < sampleInner.size(); i++)
        {
            const size_t index = sampleInner[i];
            if(mInnerExists[index])
                continue;

            vec4_t bary;
            baryComputer(mShell.mInnerShell(colon(), index), bary);

            if(min(bary) <= -1e-6)
                continue;

            double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2] + f3 * bary[3];
            double error = fabs(f + 1);
            if(maxError < error)
            {
                maxError = error;
                maxErrorPoint = PointInfo(POINT_INNER, index);
            }
            mInnerError[index] = error;
            info.points.push_back(PointInfo(POINT_INNER, index));
        }
        for(int i = 0; i < sampleOuter.size(); i++)
        {
            const size_t index = sampleOuter[i];
            if(mOuterExists[index])
                continue;

            vec4_t bary;
            baryComputer(mShell.mOuterShell(colon(), index), bary);

            if(min(bary) <= -1e-6)
                continue;

            double f = f0 * bary[0] + f1 * bary[1] + f2 * bary[2] + f3 * bary[3];
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

    void Refinement::computeAABB(const Cell &cell, double &xmax, double &xmin, double &ymax, double &ymin, double& zmax, double& zmin)
    {
        const Point& p0 = cell.vertex(0)->point();
        const Point& p1 = cell.vertex(1)->point();
        const Point& p2 = cell.vertex(2)->point();
        const Point& p3 = cell.vertex(3)->point();

        xmax = std::max(std::max(std::max(p0[0], p1[0]), p2[0]), p3[0]);
        xmin = std::min(std::min(std::min(p0[0], p1[0]), p2[0]), p3[0]);
        ymax = std::max(std::max(std::max(p0[1], p1[1]), p2[1]), p3[1]);
        ymin = std::min(std::min(std::min(p0[1], p1[1]), p2[1]), p3[1]);
        zmax = std::max(std::max(std::max(p0[2], p1[2]), p2[2]), p3[2]);
        zmin = std::min(std::min(std::min(p0[2], p1[2]), p2[2]), p3[2]);
    }

    PointType Refinement::getPointType(const VertexHandle &vh)
    {
        const PointInfo& info = vh->info();
        return info.pointType;
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
            return std::numeric_limits<double>::lowest();
        default:
            //this should not be run, otherwise there exists logic error
            throw std::logic_error("In getError(), pointType should be POINT_INNER or POINT_OUTER");
        }
    }

    void Refinement::getPointMatrix(const PointInfo &point, matrixr_t &pointMatrix)
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
            if(iterCount %100 == 0)
                std::cout << "Iteration " << iterCount << " ..." << std::endl;

            if(mSatisfyCondition1 == false)
            {
                insertPoint(mNextInsertPoint);
                updateErrors();
            }
            else
            {
                resolveBadCell();
                mBadCell = nullptr;
            }
        }

        std::cout << "Finished refinement after " << iterCount << " iterations." << std::endl;
        std::cout << "bad cell " << bad_cell << std::endl << "bad classification";
        is_outer ? std::cout << " outer " : std::cout << " inner ";
        std::cout << bad_sample_index << std::endl;
        organizeOutput();
    }

    void Refinement::organizeOutput()
    {
        //organize output data
        mOutput.vertices.resize(3, mDelaunay.number_of_vertices());
        mOutput.cells.resize(4, mDelaunay.number_of_cells());
        mOutput.vertType.reserve(mDelaunay.number_of_vertices());
        int i = 0;
        for(auto iter = mDelaunay.finite_vertices_begin(); iter != mDelaunay.finite_vertices_end(); ++iter, ++i)
        {
            mOutput.vertices(0, i) = iter->point()[0];
            mOutput.vertices(1, i) = iter->point()[1];
            mOutput.vertices(2, i) = iter->point()[2];

            PointInfo& info = iter->info();
            iter->info().indexInDelaunay = i;
            mOutput.vertType.push_back(info.pointType);
        }
        i = 0;
        for(auto iter = mDelaunay.finite_cells_begin(); iter != mDelaunay.finite_cells_end(); ++iter, ++i)
        {
            const PointInfo& info0 = iter->vertex(0)->info();
            const PointInfo& info1 = iter->vertex(1)->info();
            const PointInfo& info2 = iter->vertex(2)->info();
            const PointInfo& info3 = iter->vertex(3)->info();

            mOutput.cells(0, i) = info0.indexInDelaunay;
            mOutput.cells(1, i) = info1.indexInDelaunay;
            mOutput.cells(2, i) = info2.indexInDelaunay;
            mOutput.cells(3, i) = info3.indexInDelaunay;
        }
    }

    bool Refinement::isFinished()
    {
        //check for condition 1.
        if(mNextInsertPoint.pointType == POINT_UNKNOWN)
            return false;

        if(getError(mNextInsertPoint) > 1 - mAlpha)
        {
            mSatisfyCondition1 = false;
            return false;
        }

        mSatisfyCondition1 = true;

        //check for condition 2 and 3.
        int i = 0;
        for(auto iter = mDelaunay.finite_cells_begin(); iter != mDelaunay.finite_cells_end(); ++iter, ++i)
        {
            Cell& cell = *iter;
            const VertexHandle& vh0 = cell.vertex(0);
            const VertexHandle& vh1 = cell.vertex(1);
            const VertexHandle& vh2 = cell.vertex(2);
            const VertexHandle& vh3 = cell.vertex(3);

            if(getFValue(vh0) == getFValue(vh1) && getFValue(vh0) == getFValue(vh2) && getFValue(vh0) == getFValue(vh3))
            {
                continue;
            }

            if(cell.info().isJudged)
                continue;

            if(cell.info().notResolvable)
                continue;

            cell.info().isJudged = true;
            if(isBadCell(cell))
            {
                mBadCell = &cell;
                cell.info().isBad = true;
                bad_cell = i;
                return false;
            }
        }

        return true;
    }

    bool Refinement::isBadCell(Cell &cell)
    {
        //condition 2
        //if(computeHeight(cell) < 2 * mSampleRadius / mAlpha)
        //{
        //    return true;
        //}

        if(checkCondition3(cell) == false)
        {
            return true;
        }
        return false;
    }

    double Refinement::computeHeight(const Cell &cell)
    {
        std::vector<Point> innerVerts, outerVerts;
        for(int i = 0; i < 4; i++)
        {
            getFValue(cell.vertex(i)) < 0 ? innerVerts.push_back(cell.vertex(i)->point()) : outerVerts.push_back(cell.vertex(i)->point());
        }

        vec3_t a, b, c, d;
        if(innerVerts.size() == 1 || innerVerts.size() == 3)
        {
            if(innerVerts.size() == 1)
            {
                for(int j = 0; j < 3; j++)
                    a[j] = innerVerts[0][j];
                for(int j = 0; j < 3; j++)
                    b[j] = outerVerts[0][j];
                for(int j = 0; j < 3; j++)
                    c[j] = outerVerts[1][j];
                for(int j = 0; j < 3; j++)
                    d[j] = outerVerts[2][j];
            }
            else
            {
                for(int j = 0; j < 3; j++)
                    a[j] = outerVerts[0][j];
                for(int j = 0; j < 3; j++)
                    b[j] = innerVerts[0][j];
                for(int j = 0; j < 3; j++)
                    c[j] = innerVerts[1][j];
                for(int j = 0; j < 3; j++)
                    d[j] = innerVerts[2][j];
            }
            vec3_t bc = c - b;
            vec3_t bd = d - b;
            vec3_t h = cross(bc, bd);
            h /= norm(h);
            vec3_t ab = b - a;
            return fabs(dot(ab, h));
        }
        for(int j = 0; j < 3; j++)
            a[j] = innerVerts[0][j];
        for(int j = 0; j < 3; j++)
            b[j] = innerVerts[1][j];
        for(int j = 0; j < 3; j++)
            c[j] = outerVerts[0][j];
        for(int j = 0; j < 3; j++)
            d[j] = outerVerts[1][j];
        vec3_t n = cross(a - b, c - d);
        n /= norm(n);
        vec3_t ac = c - a;
        return fabs(dot(ac, n));
    }

    bool Refinement::checkCondition3(Cell &cell)
    {
        matrixr_t tetra(3, 4);
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                tetra(j, i) = cell.vertex(i)->point()[j];
            }
        }

        bool result = NormalChecker::check(tetra, getPointType(cell.vertex(0)), getPointType(cell.vertex(1)),
                                           getPointType(cell.vertex(2)), getPointType(cell.vertex(3)), mShell);
        return result;
    }

    bool Refinement::checkClassification(const Cell& cell, const BaryComputer &baryComputer, const matrixr_t &point, bool isOuter)
    {
        const VertexHandle& vh0 = cell.vertex(0);
        const VertexHandle& vh1 = cell.vertex(1);
        const VertexHandle& vh2 = cell.vertex(2);
        const VertexHandle& vh3 = cell.vertex(3);
        vec4_t bary;
        baryComputer(point, bary);
        double fvalue = getFValue(vh0) * bary[0] + getFValue(vh1) * bary[1] + getFValue(vh2) * bary[2] + getFValue(vh3) * bary[3];
        bool result;
        if(isOuter)
        {
            result = fvalue > 0;
        }
        else
        {
            result = fvalue < 0;
        }
        return result;
    }

    bool Refinement::resolveBadCell()
    {
        Point circumcenter = CGAL::circumcenter(mBadCell->vertex(0)->point(), mBadCell->vertex(1)->point(),
                                                mBadCell->vertex(2)->point(), mBadCell->vertex(3)->point());
        vec3_t center;
        center[0] = circumcenter[0];
        center[1] = circumcenter[1];
        center[2] = circumcenter[2];

        //compute radius of the circum circle
        vec3_t v0;
        v0[0] = mBadCell->vertex(0)->point()[0];
        v0[1] = mBadCell->vertex(0)->point()[1];
        v0[2] = mBadCell->vertex(0)->point()[2];
        double radius = norm(center - v0);

        size_t innerNearest = mShell.getInnerTree().getNearestPoint(center);
        size_t outerNearest = mShell.getOuterTree().getNearestPoint(center);
        vec3_t innerP = mShell.mInnerShell(colon(), innerNearest);
        vec3_t outerP = mShell.mOuterShell(colon(), outerNearest);        

        PointInfo info;
        vec3_t nearest_to_insert;
        if(norm(innerP - center) < norm(outerP - center))
        {
            info.pointType = POINT_INNER;
            info.index = innerNearest;
            nearest_to_insert = innerP;
        }
        else
        {
            info.pointType = POINT_OUTER;
            info.index = outerNearest;
            nearest_to_insert = outerP;
        }

        if(fabs(radius - norm(center - nearest_to_insert)) <= 1e-6)
        {
            //the nearest point is already inside the circum circle, so this bad cell is unable to resolve
            std::cerr << "[WARNING] unable to resolve bad cell" << std::endl;
            mBadCell->info().notResolvable = true;
            return false;
        }
        insertPoint(info);
        updateErrors();
        return true;
    }

    bool Refinement::insertPoint(const PointInfo &info)
    {
        switch(info.pointType)
        {
        case POINT_INNER:
            if(mInnerExists[info.index])
            {
                std::cerr << "warning : Inserting repeated inner shell point. index " << info.index << std::endl;
                return false;
            }
            mInnerExists[info.index] = true;
            mInnerError[info.index] = 0;
            break;
        case POINT_OUTER:
            if(mOuterExists[info.index])
            {
                std::cerr << "warning : Inserting repeated outer shell point. index " << info.index << std::endl;
                return false;
            }
            mOuterExists[info.index] = true;
            mOuterError[info.index] = 0;
            break;
        }

        matrixr_t point;
        getPointMatrix(info, point);
        VertexHandle vertexHandle = mDelaunay.insert(Point(point[0], point[1], point[2]));
        vertexHandle->info() = info;
        return true;
    }
}
