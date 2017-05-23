#include "EdgeCollapse.h"
#include "KernelRegion.h"
#include "SamplingQuadTree.h"
#include <wkylib/geometry.h>
#include <iostream>
#include <omp.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/io.h>

using namespace zjucad::matrix;

namespace SBV
{
    EdgeCollapse::EdgeCollapse(TriangulatedShell &triangulation, const Shell& shell, Type type, bool isHalfEdge, double sampleRadius)
        : mTriangulation(triangulation),
          mShell(shell),
          mType(type),
          mIsHalfEdge(isHalfEdge),
          mSampleRadius(sampleRadius)
    {
        mCollapseTo.reserve(triangulation.vertices.size(2));
        for(int i = 0; i < triangulation.vertices.size(2); i++)
        {
            mCollapseTo.push_back(i);
        }

        buildEdgeInfo();
    }

    void EdgeCollapse::buildEdgeInfo()
    {
        mCollapseableNeighbours.resize(mTriangulation.vertices.size(2));
        mNeighbours.resize(mTriangulation.vertices.size(2));
        mNeighbourFaces.resize(mTriangulation.vertices.size(2));
        mRelatedEdgeInfo.resize(mTriangulation.vertices.size(2));

        //build vertices neighbour info
        for(int i = 0; i < mTriangulation.triangles.size(2); i++)
        {
            size_t a = mTriangulation.triangles(0, i);
            size_t b = mTriangulation.triangles(1, i);
            size_t c = mTriangulation.triangles(2, i);

            addNeighbour(a, b);
            addNeighbour(a, c);
            addNeighbour(b, c);

            mNeighbourFaces[a].insert(i);
            mNeighbourFaces[b].insert(i);
            mNeighbourFaces[c].insert(i);

            if(mTriangulation.getSign(a) == mTriangulation.getSign(b) && mTriangulation.getSign(a) == mTriangulation.getSign(c))
            {
                continue;
            }
            if(mTriangulation.vertType[a] == POINT_BOUNDING_BOX || mTriangulation.vertType[b] == POINT_BOUNDING_BOX
                    || mTriangulation.vertType[c] == POINT_BOUNDING_BOX)
            {
                continue;
            }

            tryAddCollapseableNeighbour(a, b);
            tryAddCollapseableNeighbour(a, c);
            tryAddCollapseableNeighbour(b, c);
        }

        buildMatrices();

        //build edge infos
//#pragma omp parallel for schedule(dynamic, 1)
        for(int i = 0; i < mTriangulation.vertices.size(2); i++)
        {
            if(mTriangulation.vertType[i] == POINT_BOUNDING_BOX)
            {
                //Bounding box vertices can not be collapsed, so no need for building their edge info
                continue;
            }

            insertEdges(i);
        }
    }

    void EdgeCollapse::insertEdges(size_t vert)
    {
        //iterate over all neighbour vertices
        for(const size_t& neighbourVert : mCollapseableNeighbours[vert])
        {
            if(mIsHalfEdge)
            {
                std::shared_ptr<EdgeInfo> edgeInfo(new EdgeInfo());
                edgeInfo->firstVert = vert;
                edgeInfo->secondVert = neighbourVert;
                edgeInfo->position = mTriangulation.vertices(colon(), neighbourVert);
                edgeInfo->error = computeError(vert, edgeInfo->position) + computeError(neighbourVert, edgeInfo->position);

                mtx.lock();
                mQueue.push(edgeInfo);
                mRelatedEdgeInfo[vert].push_back(edgeInfo);
                mtx.unlock();
            }
            else
            {
                std::shared_ptr<EdgeInfo> edgeInfo(new EdgeInfo());
                edgeInfo->firstVert = vert;
                edgeInfo->secondVert = neighbourVert;
                bool found = findCollapsePos(vert, neighbourVert, edgeInfo->position, edgeInfo->error);
                if(found)
                {
                    mQueue.push(edgeInfo);
                    mRelatedEdgeInfo[vert].push_back(edgeInfo);
                }
            }
        }
    }

    void EdgeCollapse::addNeighbour(size_t firstVert, size_t secondVert)
    {
        mNeighbours[firstVert].insert(secondVert);
        mNeighbours[secondVert].insert(firstVert);
    }

    void EdgeCollapse::tryAddCollapseableNeighbour(size_t firstVert, size_t secondVert)
    {
        if(isCollapseable(firstVert, secondVert))
        {
            mCollapseableNeighbours[firstVert].insert(secondVert);
            mCollapseableNeighbours[secondVert].insert(firstVert);
        }
    }

    bool EdgeCollapse::isCollapseable(size_t firstVert, size_t secondVert)
    {
        switch(mType)
        {
        case BOUNDARY:
            return mTriangulation.vertType[firstVert] == mTriangulation.vertType[secondVert]
                    && mTriangulation.vertType[firstVert] != POINT_BOUNDING_BOX
                    && mTriangulation.vertType[firstVert] != POINT_ZERO;
        case ZERO_SET:
            return mTriangulation.vertType[firstVert] == mTriangulation.vertType[secondVert]
                    && mTriangulation.vertType[firstVert] == POINT_ZERO;
        default:
            throw std::logic_error("not implemented.");
        }
    }

    void EdgeCollapse::buildMatrices()
    {
        //initialize error matrices to zero
        for(int i = 0; i < mTriangulation.vertices.size(2); i++)
        {
            mQ.push_back(Eigen::Matrix3d::Zero());
        }

        for(int i = 0; i < mTriangulation.vertices.size(2); i++)
        {
            if(mTriangulation.vertType[i] == POINT_BOUNDING_BOX)
            {
                //Bounding box vertices can not be collapsed, so no need for computing their error matrices
                continue;
            }

            const matrixr_t&  posI = mTriangulation.vertices(colon(), i);
            auto& neighbour = mCollapseableNeighbours[i];
            for(const size_t& neighbourVert : neighbour)
            {
                Eigen::Vector3d p;

                matrixr_t a = posI - mTriangulation.vertices(colon(), neighbourVert);
                a = a / norm(a);
                p[0] = a[1];
                p[1] = -a[0];
                p[2] = posI[1] * a[0] - posI[0] * a[1];

                mQ[i] += p * p.transpose();
            }
        }
    }

    double EdgeCollapse::computeError(size_t vert, const matrixr_t &point)
    {
        Eigen::Vector3d p;
        p[0] = point[0];
        p[1] = point[1];
        p[2] = 1;
        return p.transpose() * mQ[vert] * p;
    }

    void EdgeCollapse::collapse()
    {
        int numCollapsed = 0;
        int i = 0;
        while(!mQueue.empty())
        {
            i++;
            std::shared_ptr<EdgeInfo> edgeInfo = mQueue.top();
            mQueue.pop();
            if(edgeInfo->isUpdated)
            {
                //this edge info is out-dated
                continue;
            }
            if((mCollapseTo[edgeInfo->firstVert] != edgeInfo->firstVert)
                    || (mCollapseTo[edgeInfo->secondVert] != edgeInfo->secondVert))
            {
                //the vertices of this edge has been collapsed
                continue;
            }

            if(!isValidCollapse(edgeInfo->firstVert, edgeInfo->secondVert, mTriangulation.vertices(colon(), edgeInfo->secondVert)))
            {
                continue;
            }

            collapseEdge(edgeInfo->firstVert, edgeInfo->secondVert, edgeInfo->position);
            numCollapsed++;
            std::cout << "Collapse " << numCollapsed << " : from " << edgeInfo->firstVert << " to " << edgeInfo->secondVert << std::endl;
        }

        organizeOutput();
    }

    void EdgeCollapse::organizeOutput()
    {
        std::vector<int> finalIndex;   //recording the final vert index of the original vertices
        finalIndex.reserve(mTriangulation.vertices.size(2));
        for(int i = 0; i < mTriangulation.vertices.size(2); i++)
        {
            finalIndex.push_back(-1);
        }

        //compute output vertices count
        size_t newVertCount = 0;
        for(int i = 0; i < mTriangulation.vertices.size(2); i++)
        {
            if(mCollapseTo[i] == i)
            {
                newVertCount++;
            }
        }

        matrixr_t newVertices(2, newVertCount);

        //organize new vertices
        for(int i = 0, j = 0; i < mTriangulation.vertices.size(2); i++)
        {
            size_t collapsedVert = getCollapsedVert(i);
            if(finalIndex[i] >= 0)
            {
                //the vert has already been dispatched, so continue
                continue;
            }
            if(finalIndex[collapsedVert] >= 0)
            {
                //the final vertex has been added
                finalIndex[i] = finalIndex[collapsedVert];
                continue;
            }
            finalIndex[i] = j;
            if(collapsedVert == i)
            {
                //the vertex is not collapsed, so no change needed
                newVertices(colon(), j++) = mTriangulation.vertices(colon(), i);
            }
            else
            {
                finalIndex[collapsedVert] = j;
                newVertices(colon(), j++) = mTriangulation.vertices(colon(), collapsedVert);
            }
        }

        //update F value
        for(int i = 0; i < mTriangulation.vertices.size(2); i++)
        {
            mTriangulation.vertType[finalIndex[i]] = mTriangulation.vertType[i];
        }
        mTriangulation.vertType.erase(mTriangulation.vertType.begin() + newVertCount, mTriangulation.vertType.end());

        //organize new triangles, remove duplicate
        std::vector<matrixs_t> newTriangleVector;
        for(int i = 0; i < mTriangulation.triangles.size(2); i++)
        {
            bool isDuplicated = false;
            std::set<size_t> triangle;

            size_t a = finalIndex[mTriangulation.triangles(0, i)];
            size_t b = finalIndex[mTriangulation.triangles(1, i)];
            size_t c = finalIndex[mTriangulation.triangles(2, i)];

            if(a == b || a == c || b == c)
            {
                //the face is collapsed
                continue;
            }

            triangle.insert(a);
            triangle.insert(b);
            triangle.insert(c);
            for(int j = 0; j < newTriangleVector.size(); j++)
            {
                std::set<size_t> oldTriangle;
                oldTriangle.insert(newTriangleVector[j][0]);
                oldTriangle.insert(newTriangleVector[j][1]);
                oldTriangle.insert(newTriangleVector[j][2]);
                if(triangle == oldTriangle)
                {
                    isDuplicated = true;
                    break;
                }
            }
            if(isDuplicated == false)
            {
                matrixs_t tri(3, 1);
                tri[0] = a;
                tri[1] = b;
                tri[2] = c;
                newTriangleVector.push_back(tri);
            }
        }

        matrixs_t newTriangles(3, newTriangleVector.size());
        for(int i = 0; i < newTriangleVector.size(); i++)
        {
            newTriangles(colon(), i) = newTriangleVector[i];
        }

        mTriangulation.vertices = newVertices;
        mTriangulation.triangles = newTriangles;
    }

    void EdgeCollapse::collapseEdge(size_t firstVert, size_t secondVert, const matrixr_t &collapseTo)
    {
        mTriangulation.vertices(colon(), secondVert) = collapseTo;
        updateEdgeInfo(firstVert, secondVert);
    }

    void EdgeCollapse::mergeNeighbours(size_t vertCollapsed, size_t vertCollapsedTo)
    {
        //merge collapseable neighbours
        for(size_t vert : mCollapseableNeighbours[vertCollapsed])
        {
            if(vert == vertCollapsedTo)
            {
                mCollapseableNeighbours[vert].erase(mCollapseableNeighbours[vert].find(vertCollapsed));
                continue;
            }
            else
            {
                mCollapseableNeighbours[vert].erase(mCollapseableNeighbours[vert].find(vertCollapsed));
                mCollapseableNeighbours[vert].insert(vertCollapsedTo);
            }
            mCollapseableNeighbours[vertCollapsedTo].insert(vert);
        }

        //merge neighbours
        for(size_t vert : mNeighbours[vertCollapsed])
        {
            if(vert == vertCollapsedTo)
            {
                mNeighbours[vert].erase(mNeighbours[vert].find(vertCollapsed));
                continue;
            }
            else
            {
                mNeighbours[vert].erase(mNeighbours[vert].find(vertCollapsed));
                mNeighbours[vert].insert(vertCollapsedTo);
            }
            mNeighbours[vertCollapsedTo].insert(vert);
        }

        //merge neighbour faces
        std::set<size_t>& nfCollapsedTo = mNeighbourFaces[vertCollapsedTo];
        for(size_t face : mNeighbourFaces[vertCollapsed])
        {
            auto iter = nfCollapsedTo.find(face);
            if(iter != nfCollapsedTo.end())
            {
                nfCollapsedTo.erase(iter);
            }
            else
            {
                nfCollapsedTo.insert(face);
            }
        }
    }

    void EdgeCollapse::updateEdgeInfo(size_t vertCollapsed, size_t vertCollapsedTo)
    {
        mCollapseTo[vertCollapsed] = vertCollapsedTo;
        mQ[vertCollapsedTo] += mQ[vertCollapsed];

        //mark related edge info as out-dated
        for(auto& edgeInfoOutdated : mRelatedEdgeInfo[vertCollapsedTo])
        {
            edgeInfoOutdated->isUpdated = true;
        }

        //now the edge infos has been marked, so clear the vector
        mRelatedEdgeInfo[vertCollapsedTo].clear();

        mergeNeighbours(vertCollapsed, vertCollapsedTo);

        //insert new edge infos
        insertEdges(vertCollapsedTo);
    }

    size_t EdgeCollapse::getCollapsedVert(size_t vert)
    {
        while(mCollapseTo[vert] != vert)
        {
            vert = mCollapseTo[vert];
        }
        return vert;
    }

    bool EdgeCollapse::isValidCollapse(size_t firstVert, size_t secondVert, const matrixr_t &collapseTo)
    {
        if(testLinkCondition(firstVert, secondVert) == false)
        {
            return false;
        }

        matrixs_t lines;
        std::set<size_t> innerSample;
        std::set<size_t> outerSample;
        buildOneRingArea(firstVert, secondVert, lines, innerSample, outerSample);

        KernelRegion kernel(mTriangulation.vertices, lines, mShell, innerSample, outerSample,
                            mTriangulation, mTriangulation.vertType[firstVert]);
        if(kernel.contains(collapseTo) == false)
        {
            return false;
        }

        return true;
    }

    bool EdgeCollapse::testLinkCondition(size_t firstVert, size_t secondVert)
    {
        //we mark A as firstVert, B as secondVert in this function
        const std::set<size_t>& neighbourA = mNeighbours[firstVert];
        const std::set<size_t>& neighbourB = mNeighbours[secondVert];

        std::set<size_t> linkAB;

        //build linkAB, which is common neighbour verts shared by A and B
        for(auto& neighbour : neighbourA)
        {
            size_t collapsedNeighbour = getCollapsedVert(neighbour);
            for(auto& neighbour2 : neighbourB)
            {
                size_t collapsedNeighbour2 = getCollapsedVert(neighbour2);
                if(collapsedNeighbour == collapsedNeighbour2)
                {
                    linkAB.insert(neighbour);
                    break;
                }
            }
        }

        //test whether each vert in linkAB can form a face with A and B
        for(auto& neighbour : linkAB)
        {
            bool isFace = false;
            for(auto& face : mNeighbourFaces[firstVert])
            {
                matrixs_t f = mTriangulation.triangles(colon(), face);
                size_t a = getCollapsedVert(f[0]);
                size_t b = getCollapsedVert(f[1]);
                size_t c = getCollapsedVert(f[2]);

                std::set<size_t> triangle1;
                std::set<size_t> triangle2;
                triangle1.insert(a);
                triangle1.insert(b);
                triangle1.insert(c);
                triangle2.insert(firstVert);
                triangle2.insert(secondVert);
                triangle2.insert(neighbour);

                if(triangle1 == triangle2)
                {
                    isFace = true;
                    break;
                }
            }
            if(isFace == false)
            {
                //If there is a vert in linkAB that cannot form a face with A and B, then link condition test fails
                return false;
            }
        }

        //all verts in linkAB can form a face with A and B, so link condition test passes
        return true;
    }

    void EdgeCollapse::buildOneRingArea(size_t firstVert, size_t secondVert, matrixs_t& lines,
                                            std::set<size_t>& innerSample, std::set<size_t>& outerSample)
    {
        std::vector<std::pair<size_t, size_t> > boundaryEdges;
        findBoundaryEdge(firstVert, secondVert, boundaryEdges);
        findBoundaryEdge(secondVert, firstVert, boundaryEdges);
        findShellSamples(firstVert, innerSample, outerSample);
        findShellSamples(secondVert, innerSample, outerSample);

        lines.resize(2, boundaryEdges.size());
        int i = 0;
        for(std::pair<size_t, size_t>& edge : boundaryEdges)
        {
            lines(0, i) = edge.first;
            lines(1, i) = edge.second;
            i++;
        }
    }

    void EdgeCollapse::findBoundaryEdge(size_t firstVert, size_t secondVert, std::vector<std::pair<size_t, size_t> >& boundaryEdges)
    {
        for(size_t face : mNeighbourFaces[firstVert])
        {
            size_t a = getCollapsedVert(mTriangulation.triangles(0, face));
            size_t b = getCollapsedVert(mTriangulation.triangles(1, face));
            size_t c = getCollapsedVert(mTriangulation.triangles(2, face));

            if(a == b || a == c || b == c)
            {
                //the face is collapsed
                continue;
            }

            if(a == secondVert || b == secondVert || c == secondVert)
            {
                continue;
            }

            if(a == firstVert)
            {
                boundaryEdges.push_back(std::make_pair(b, c));
            }
            else if(b == firstVert)
            {
                boundaryEdges.push_back(std::make_pair(a, c));
            }
            else
            {
                boundaryEdges.push_back(std::make_pair(a, b));
            }
        }
    }

    void EdgeCollapse::findShellSamples(size_t vert, std::set<size_t> &innerSample, std::set<size_t> &outerSample)
    {
        double xmin = std::numeric_limits<double>::max();
        double xmax = std::numeric_limits<double>::min();
        double ymin = std::numeric_limits<double>::max();
        double ymax = std::numeric_limits<double>::min();

        // AABB
        for(size_t face : mNeighbourFaces[vert])
        {
            const size_t a = getCollapsedVert(mTriangulation.triangles(0, face));
            const size_t b = getCollapsedVert(mTriangulation.triangles(1, face));
            const size_t c = getCollapsedVert(mTriangulation.triangles(2, face));

            if(a == b || a == c || b == c)
            {
                //the face is collapsed
                continue;
            }

            for(int i = 0; i < 3; i++)
            {
                const matrixr_t& vert = mTriangulation.vertices(colon(), getCollapsedVert(mTriangulation.triangles(i, face)));
                if(vert[0] > xmax)
                {
                    xmax = vert[0];
                }
                if(vert[0] < xmin)
                {
                    xmin = vert[0];
                }
                if(vert[1] > ymax)
                {
                    ymax = vert[1];
                }
                if(vert[1] < ymin)
                {
                    ymin = vert[1];
                }
            }
        }

        matrixs_t sampleInner;
        mShell.getInnerTree().getPointsInRange(xmin, xmax, ymin, ymax, sampleInner);
        matrixs_t sampleOuter;
        mShell.getOuterTree().getPointsInRange(xmin, xmax, ymin, ymax, sampleOuter);

        for(size_t face : mNeighbourFaces[vert])
        {
            size_t a = getCollapsedVert(mTriangulation.triangles(0, face));
            size_t b = getCollapsedVert(mTriangulation.triangles(1, face));
            size_t c = getCollapsedVert(mTriangulation.triangles(2, face));

            if(a == b || a == c || b == c)
            {
                //the face is collapsed, or the sample point is the vert of triangle
                continue;
            }

            matrixr_t invA = ones<double>(3, 3);
            invA(colon(0, 1), 0) = mTriangulation.vertices(colon(), a);
            invA(colon(0, 1), 1) = mTriangulation.vertices(colon(), b);
            invA(colon(0, 1), 2) = mTriangulation.vertices(colon(), c);

            if(inv(invA)) {
                std::cerr << "warning: degenerated triangle." << std::endl;
                continue;
            }
            matrixr_t barysInner = invA(colon(), colon(0, 1)) * mShell.mInnerShell(colon(), sampleInner) +
                    invA(colon(), 2) * ones<double>(1, sampleInner.size());
            for(int i = 0; i < sampleInner.size(); i++)
            {
                if(min(barysInner(colon(), i)) >= 0)
                {
                    innerSample.insert(sampleInner[i]);
                }
            }
            matrixr_t barysOuter = invA(colon(), colon(0, 1)) * mShell.mOuterShell(colon(), sampleOuter) +
                    invA(colon(), 2) * ones<double>(1, sampleOuter.size());
            for(int i = 0; i < sampleOuter.size(); i++)
            {
                if(min(barysOuter(colon(), i)) >= 0)
                {
                    outerSample.insert(sampleOuter[i]);
                }
            }

            /*for(int i = 0; i < sampleInner.size(); i++)
            {
                const matrixr_t b = invA(colon(), colon(0, 1))* mShell.mInnerShell(colon(), sampleInner[i])
                        + invA(colon(), 2);
                if(min(b) >= 0)
                {
                    innerSample.insert(sampleInner[i]);
                }
            }
            for(int i = 0; i < sampleOuter.size(); i++)
            {
                const matrixr_t b = invA(colon(), colon(0, 1))* mShell.mOuterShell(colon(), sampleOuter[i])
                        + invA(colon(), 2);
                if(min(b) >= 0)
                {
                    outerSample.insert(sampleOuter[i]);
                }
            }*/

            /*matrixr_t triangle = mTriangulation.vertices(colon(), mTriangulation.triangles(colon(), face));
            triangle(colon(), 0) = mTriangulation.vertices(colon(), a);
            triangle(colon(), 1) = mTriangulation.vertices(colon(), b);
            triangle(colon(), 2) = mTriangulation.vertices(colon(), c);

            for(int i = 0; i < sampleInner.size(); i++)
            {
                const matrixr_t& point = mShell.mInnerShell(colon(), sampleInner[i]);

                if(fabs(max(triangle(colon(), 0) - point)) < 1e-6
                        || fabs(max(triangle(colon(), 1) - point)) < 1e-6
                        || fabs(max(triangle(colon(), 2) - point)) < 1e-6 )
                {
                    continue;
                }
                matrixr_t bary;
                if(WKYLIB::barycentric_2D(point, triangle, bary))
                {
                    innerSample.insert(sampleInner[i]);
                }
            }

            for(int i = 0; i < sampleOuter.size(); i++)
            {
                const matrixr_t& point = mShell.mOuterShell(colon(), sampleOuter[i]);

                if(fabs(max(triangle(colon(), 0) - point)) < 1e-6\
                        || fabs(max(triangle(colon(), 1) - point)) < 1e-6
                        || fabs(max(triangle(colon(), 2) - point)) < 1e-6 )
                {
                    //continue;
                }
                matrixr_t bary;
                if(WKYLIB::barycentric_2D(point, triangle, bary))
                {
                    outerSample.insert(sampleOuter[i]);
               }
            }*/
        }
    }

    bool EdgeCollapse::findCollapsePos(size_t vert, size_t vertCollapseTo, matrixr_t &position, double& out_error)
    {
        matrixs_t lines;
        std::set<size_t> innerSample;
        std::set<size_t> outerSample;
        bool found = false;
        buildOneRingArea(vert, vertCollapseTo, lines, innerSample, outerSample);

        KernelRegion kernel(mTriangulation.vertices, lines, mShell, innerSample, outerSample,
                                 mTriangulation, mTriangulation.vertType[vert]);
        out_error = std::numeric_limits<double>::max();

        if(mType == BOUNDARY)
        {
            if(mTriangulation.vertType[vert] == POINT_INNER)
            {
                for(size_t sample : innerSample)
                {
                    const matrixr_t samplePoint = mShell.mInnerShell(colon(), sample);
                    if(kernel.contains(samplePoint))
                    {
                        double error = computeError(vert, samplePoint) + computeError(vertCollapseTo, samplePoint);
                        if(error < out_error)
                        {
                            found = true;
                            out_error = error;
                            position = samplePoint;
                        }
                    }
                }
            }
            else if(mTriangulation.vertType[vert] == POINT_OUTER)
            {
                for(size_t sample : outerSample)
                {
                    const matrixr_t samplePoint = mShell.mOuterShell(colon(), sample);
                    if(kernel.contains(samplePoint))
                    {
                        double error = computeError(vert, samplePoint) + computeError(vertCollapseTo, samplePoint);
                        if(error < out_error)
                        {
                            found = true;
                            out_error = error;
                            position = samplePoint;
                        }
                    }
                }
            }
            else
            {
                throw std::logic_error("this should not be run");
            }
        }
        else if(mType == ZERO_SET)
        {
            //find AABB of the one ring area
            double xmin = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::min();
            double ymin = std::numeric_limits<double>::max();
            double ymax = std::numeric_limits<double>::min();

            for(int i = 0; i < lines.size(2); i++)
            {
                const matrixr_t& a = mTriangulation.vertices(colon(), lines(0, i));
                const matrixr_t& b = mTriangulation.vertices(colon(), lines(1, i));
                xmax = a[0] > xmax ? a[0] : xmax;
                xmin = a[0] < xmin ? a[0] : xmin;
                ymax = a[1] > ymax ? a[1] : ymax;
                ymin = a[1] < ymin ? a[1] : ymin;
                xmax = b[0] > xmax ? b[0] : xmax;
                xmin = b[0] < xmin ? b[0] : xmin;
                ymax = b[1] > ymax ? b[1] : ymax;
                ymin = b[1] < ymin ? b[1] : ymin;
            }

            SamplingQuadTree tree(kernel, xmax, xmin, ymax, ymin, mSampleRadius);
            auto& samples = tree.getSamples();
            for(int i = 0; i < samples.size(); i++)
            {
                auto& point = samples[i];
                double error = computeError(vert, point) + computeError(vertCollapseTo, point);
                if(error < out_error)
                {
                    found = true;
                    out_error = error;
                    position = point;
                }
            }
        }
        return found;
    }
}
