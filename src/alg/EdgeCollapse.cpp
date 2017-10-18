#include "EdgeCollapse.h"
#include "BaryComputer.h"
#include "KernelRegion.h"
#include "SamplingQuadTree.h"
#include "NormalChecker.h"
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
        mNeighbourCells.resize(mTriangulation.vertices.size(2));
        mSameShellNeighbourFaces.resize(mTriangulation.vertices.size(2));
        mRelatedEdgeInfo.resize(mTriangulation.vertices.size(2));

        //build vertices neighbour info
        for(int i = 0; i < mTriangulation.cells.size(2); i++)
        {
            size_t a = mTriangulation.cells(0, i);
            size_t b = mTriangulation.cells(1, i);
            size_t c = mTriangulation.cells(2, i);
            size_t d = mTriangulation.cells(3, i);

            addNeighbour(a, b);
            addNeighbour(a, c);
            addNeighbour(a, d);
            addNeighbour(b, c);
            addNeighbour(b, d);
            addNeighbour(c, d);

            mNeighbourCells[a].insert(i);
            mNeighbourCells[b].insert(i);
            mNeighbourCells[c].insert(i);
            mNeighbourCells[d].insert(i);

            if(mTriangulation.getSign(a) == mTriangulation.getSign(b) && mTriangulation.getSign(a) == mTriangulation.getSign(c)
                    && mTriangulation.getSign(a) == mTriangulation.getSign(d))
            {
                continue;
            }
            if(mTriangulation.vertType[a] == POINT_BOUNDING_BOX || mTriangulation.vertType[b] == POINT_BOUNDING_BOX
                    || mTriangulation.vertType[c] == POINT_BOUNDING_BOX || mTriangulation.vertType[d] == POINT_BOUNDING_BOX)
            {
                continue;
            }

            tryAddCollapseableNeighbour(a, b);
            tryAddCollapseableNeighbour(a, c);
            tryAddCollapseableNeighbour(a, d);
            tryAddCollapseableNeighbour(b, c);
            tryAddCollapseableNeighbour(b, d);
            tryAddCollapseableNeighbour(c, d);

            tryAddCollapseableNeighbourFace(a, b, c);
            tryAddCollapseableNeighbourFace(a, b, d);
            tryAddCollapseableNeighbourFace(a, c, d);
            tryAddCollapseableNeighbourFace(b, c, d);
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
                    mtx.lock();
                    mQueue.push(edgeInfo);
                    mRelatedEdgeInfo[vert].push_back(edgeInfo);
                    mtx.unlock();
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

    void EdgeCollapse::tryAddCollapseableNeighbourFace(size_t firstVert, size_t secondVert, size_t thirdVert)
    {
        if(isSameShellFace(firstVert, secondVert, thirdVert))
        {
            mSameShellNeighbourFaces[firstVert].insert(std::make_pair(secondVert, thirdVert));
            mSameShellNeighbourFaces[secondVert].insert(std::make_pair(firstVert, thirdVert));
            mSameShellNeighbourFaces[thirdVert].insert(std::make_pair(firstVert, secondVert));
        }
    }

    bool EdgeCollapse::isSameShellFace(size_t firstVert, size_t secondVert, size_t thirdVert)
    {
        if(mTriangulation.vertType[firstVert] == mTriangulation.vertType[secondVert] &&
                mTriangulation.vertType[firstVert] == mTriangulation.vertType[thirdVert])
        {
            return (mTriangulation.vertType[firstVert] == POINT_INNER) || (mTriangulation.vertType[secondVert] == POINT_OUTER);
        }
        return false;
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
            mQ.push_back(zeros<double>(4, 4));
        }

        for(int i = 0; i < mTriangulation.vertices.size(2); i++)
        {
            if(mTriangulation.vertType[i] == POINT_BOUNDING_BOX)
            {
                //Bounding box vertices can not be collapsed, so no need for computing their error matrices
                continue;
            }

            const vec3_t a = mTriangulation.vertices(colon(), i);
            auto& neighbourFaces = mSameShellNeighbourFaces[i];
            for(const std::pair<size_t, size_t>& neighbourFace : neighbourFaces)
            {
                vec4_t p;

                vec3_t b = mTriangulation.vertices(colon(), neighbourFace.first);
                vec3_t c = mTriangulation.vertices(colon(), neighbourFace.second);
                vec3_t ab = b - a;
                vec3_t ac = c - a;
                vec3_t n = cross(ab, ac);

                n = n / norm(n);

                p[0] = n[0];
                p[1] = n[1];
                p[2] = n[2];
                p[3] = -dot(n, a);

                mQ[i] += p * trans(p);
            }
        }
    }

    double EdgeCollapse::computeError(size_t vert, const matrixr_t &point)
    {
        matrixr_t v(4, 1);
        v[0] = point[0];
        v[1] = point[1];
        v[2] = point[2];
        v[3] = 1;
        matrixr_t temp = trans(v) * mQ[vert];
        return dot(temp, v);
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

            if(!isValidCollapse(edgeInfo->firstVert, edgeInfo->secondVert, edgeInfo->position))
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

        matrixr_t newVertices(3, newVertCount);

        //organize new vertices
        for(int i = 0, j = 0; i < mTriangulation.vertices.size(2); i++)
        {
            size_t collapsedVert = getCollapsedVert(i);
            if(finalIndex[i] >= 0)
            {
                //the vert has already been processed, so continue
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
        std::vector<matrixs_t> newCellVector;
        for(int i = 0; i < mTriangulation.cells.size(2); i++)
        {
            bool isDuplicated = false;
            std::set<size_t> cell;

            size_t a = finalIndex[mTriangulation.cells(0, i)];
            size_t b = finalIndex[mTriangulation.cells(1, i)];
            size_t c = finalIndex[mTriangulation.cells(2, i)];
            size_t d = finalIndex[mTriangulation.cells(3, i)];

            if(a == b || a == c || a== d || b == c || b == d || c == d)
            {
                //the cell is collapsed
                continue;
            }

            cell.insert(a);
            cell.insert(b);
            cell.insert(c);
            cell.insert(d);
            for(int j = 0; j < newCellVector.size(); j++)
            {
                std::set<size_t> oldCell;
                oldCell.insert(newCellVector[j][0]);
                oldCell.insert(newCellVector[j][1]);
                oldCell.insert(newCellVector[j][2]);
                oldCell.insert(newCellVector[j][3]);
                if(cell == oldCell)
                {
                    isDuplicated = true;
                    break;
                }
            }
            if(isDuplicated == false)
            {
                matrixs_t cel(4, 1);
                cel[0] = a;
                cel[1] = b;
                cel[2] = c;
                cel[3] = d;
                newCellVector.push_back(cel);
            }
        }

        matrixs_t newCells(4, newCellVector.size());
        for(int i = 0; i < newCellVector.size(); i++)
        {
            newCells(colon(), i) = newCellVector[i];
        }

        mTriangulation.vertices = newVertices;
        mTriangulation.cells = newCells;
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

        //merge neighbour cells
        std::set<size_t>& nfCollapsedTo = mNeighbourCells[vertCollapsedTo];
        for(size_t face : mNeighbourCells[vertCollapsed])
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

        //merge same shell neighbour faces
        /*std::set<std::pair<size_t, size_t>>& sameFacesCollapsedTo = mSameShellNeighbourFaces[vertCollapsedTo];
        for(auto iter = sameFacesCollapsedTo.begin(); iter != sameFacesCollapsedTo.end(); ++iter)
        {
            const std::pair<size_t, size_t>& face = *iter;
            if(getCollapsedVert(face.first) == getCollapsedVert(vertCollapsed)
                    || getCollapsedVert(face.second) == getCollapsedVert(vertCollapsed))
            {
                iter = sameFacesCollapsedTo.erase(iter);
                --iter;
            }
        }
        for(const std::pair<size_t, size_t>& face : mSameShellNeighbourFaces[vertCollapsed])
        {
            if(getCollapsedVert(face.first) != getCollapsedVert(vertCollapsedTo)
                    && getCollapsedVert(face.second) != getCollapsedVert(vertCollapsedTo))
            {
                sameFacesCollapsedTo.insert(face);
            }
        }*/
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

    bool EdgeCollapse::isValidCollapse(size_t firstVert, size_t secondVert, const vec3_t &collapseTo)
    {
        if(testLinkCondition(firstVert, secondVert) == false)
        {
            return false;
        }

        matrixs_t faces;
        std::set<size_t> innerSample;
        std::set<size_t> outerSample;
        buildOneRingArea(firstVert, secondVert, faces, innerSample, outerSample);

        KernelRegion kernel(mTriangulation.vertices, faces, mTriangulation.vertices(colon(), firstVert),
                            mShell, innerSample, outerSample,
                            mTriangulation, mTriangulation.vertType[firstVert]);
        if(kernel.contains(collapseTo) == false)
        {
            return false;
        }

        if(checkNormal(faces, collapseTo, mTriangulation.vertType[firstVert]) == false)
            return false;

        return true;
    }

    bool EdgeCollapse::testLinkCondition(size_t firstVert, size_t secondVert)
    {
        //we mark A as firstVert, B as secondVert in this function
        const std::set<size_t>& neighbourA = mNeighbourCells[firstVert];
        const std::set<size_t>& neighbourB = mNeighbourCells[secondVert];

        //build linkAB, which is common neighbour edges shared by A and B
        //"neighbour edge" is defined as below:
        //    B
        //   /|\
        //  A-|-D
        //   \|/
        //    C
        //for this tetrahedron, the "neighbour edge" of A is BC, BD and CD.
        using Edge = std::set<size_t>;

        std::set<Edge> neighbourEdgesA;
        for(auto& cell : neighbourA)
        {
            std::vector<size_t> offVerts;
            for(int i = 0; i < 4; i++)
            {
                size_t vert = getCollapsedVert(mTriangulation.cells(i, cell));
                if(vert != firstVert)
                    offVerts.push_back(vert);
            }
            Edge edge1, edge2, edge3;
            edge1.insert(offVerts[0]);
            edge1.insert(offVerts[1]);
            edge2.insert(offVerts[0]);
            edge2.insert(offVerts[2]);
            edge3.insert(offVerts[1]);
            edge3.insert(offVerts[2]);
            neighbourEdgesA.insert(edge1);
            neighbourEdgesA.insert(edge2);
            neighbourEdgesA.insert(edge3);
        }

        std::set<Edge> neighbourEdgesB;
        for(auto& cell : neighbourB)
        {
            std::vector<size_t> offVerts;
            for(int i = 0; i < 4; i++)
            {
                size_t vert = getCollapsedVert(mTriangulation.cells(i, cell));
                if(vert != secondVert)
                    offVerts.push_back(vert);
            }
            Edge edge1, edge2, edge3;
            edge1.insert(offVerts[0]);
            edge1.insert(offVerts[1]);
            edge2.insert(offVerts[0]);
            edge2.insert(offVerts[2]);
            edge3.insert(offVerts[1]);
            edge3.insert(offVerts[2]);
            neighbourEdgesB.insert(edge1);
            neighbourEdgesB.insert(edge2);
            neighbourEdgesB.insert(edge3);
        }

        //build linkAB, which is common neighbour edges of A and B.
        std::set<Edge> linkAB;
        for(auto& edge : neighbourEdgesA)
        {
            if(neighbourEdgesB.find(edge) != neighbourEdgesB.end())
            {
                linkAB.insert(edge);
            }
        }

        //test whether each edge in linkAB can form a tetrahedron with A and B
        for(auto& neighbourEdge : linkAB)
        {
            bool isTetra = false;
            for(auto& cell : mNeighbourCells[firstVert])
            {
                matrixs_t f = mTriangulation.cells(colon(), cell);
                size_t a = getCollapsedVert(f[0]);
                size_t b = getCollapsedVert(f[1]);
                size_t c = getCollapsedVert(f[2]);
                size_t d = getCollapsedVert(f[3]);

                std::set<size_t> tetra1;
                std::set<size_t> tetra2;
                tetra1.insert(a);
                tetra1.insert(b);
                tetra1.insert(c);
                tetra1.insert(d);
                tetra2.insert(firstVert);
                tetra2.insert(secondVert);
                for(const size_t& neighbourEdgeVert : neighbourEdge)
                {
                    tetra2.insert(neighbourEdgeVert);
                }
                if(tetra1 == tetra2)
                {
                    isTetra = true;
                    break;
                }
            }
            if(isTetra == false)
            {
                //If there is an edge in linkAB that cannot form a tetrahedron with A and B, then link condition test fails
                return false;
            }
        }

        //all edges in linkAB can form a tetrahedron with A and B, so link condition test passes
        return true;
    }

    void EdgeCollapse::buildOneRingArea(size_t firstVert, size_t secondVert, matrixs_t& faces,
                                            std::set<size_t>& innerSample, std::set<size_t>& outerSample)
    {
        std::vector<std::tuple<size_t, size_t, size_t> > boundaryFaces;
        findBoundaryFace(firstVert, secondVert, boundaryFaces);
        findBoundaryFace(secondVert, firstVert, boundaryFaces);
        findShellSamples(firstVert, innerSample, outerSample);
        findShellSamples(secondVert, innerSample, outerSample);

        faces.resize(3, boundaryFaces.size());
        int i = 0;
        for(std::tuple<size_t, size_t, size_t>& face : boundaryFaces)
        {
            faces(0, i) = std::get<0>(face);
            faces(1, i) = std::get<1>(face);
            faces(2, i) = std::get<2>(face);
            i++;
        }
    }

    void EdgeCollapse::findBoundaryFace(size_t firstVert, size_t secondVert, std::vector<std::tuple<size_t, size_t, size_t> >& boundaryFaces)
    {
        for(size_t cell : mNeighbourCells[firstVert])
        {
            size_t a = getCollapsedVert(mTriangulation.cells(0, cell));
            size_t b = getCollapsedVert(mTriangulation.cells(1, cell));
            size_t c = getCollapsedVert(mTriangulation.cells(2, cell));
            size_t d = getCollapsedVert(mTriangulation.cells(3, cell));

            if(a == b || a == c || a == d || b == c || b == d || c == d)
            {
                //the cell is collapsed
                continue;
            }

            if(a == secondVert || b == secondVert || c == secondVert || d == secondVert)
            {
                continue;
            }

            if(a == firstVert)
            {
                boundaryFaces.push_back(std::make_tuple(b, c, d));
            }
            else if(b == firstVert)
            {
                boundaryFaces.push_back(std::make_tuple(a, c, d));
            }
            else if(c == firstVert)
            {
                boundaryFaces.push_back(std::make_tuple(a, b, d));
            }
            else
            {
                boundaryFaces.push_back(std::make_tuple(a, b, c));
            }
        }
    }

    void EdgeCollapse::findShellSamples(size_t vert, std::set<size_t> &innerSample, std::set<size_t> &outerSample)
    {
        double xmin = std::numeric_limits<double>::max();
        double xmax = std::numeric_limits<double>::lowest();
        double ymin = std::numeric_limits<double>::max();
        double ymax = std::numeric_limits<double>::lowest();
        double zmin = std::numeric_limits<double>::max();
        double zmax = std::numeric_limits<double>::lowest();

        // AABB
        for(size_t cell : mNeighbourCells[vert])
        {
            size_t a = getCollapsedVert(mTriangulation.cells(0, cell));
            size_t b = getCollapsedVert(mTriangulation.cells(1, cell));
            size_t c = getCollapsedVert(mTriangulation.cells(2, cell));
            size_t d = getCollapsedVert(mTriangulation.cells(3, cell));

            if(a == b || a == c || a == d || b == c || b == d || c == d)
            {
                //the cell is collapsed
                continue;
            }

            //find AABB
            for(int i = 0; i < 4; i++)
            {
                const vec3_t vert = mTriangulation.vertices(colon(), getCollapsedVert(mTriangulation.cells(i, cell)));
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
                if(vert[2] > zmax)
                {
                    zmax = vert[2];
                }
                if(vert[2] < zmin)
                {
                    zmin = vert[2];
                }
            }
        }

        matrixs_t sampleInner;
        mShell.getInnerTree().getPointsInRange(xmin, xmax, ymin, ymax, zmin, zmax, sampleInner);
        matrixs_t sampleOuter;
        mShell.getOuterTree().getPointsInRange(xmin, xmax, ymin, ymax, zmin, zmax, sampleOuter);

        for(size_t cell : mNeighbourCells[vert])
        {
            size_t a = getCollapsedVert(mTriangulation.cells(0, cell));
            size_t b = getCollapsedVert(mTriangulation.cells(1, cell));
            size_t c = getCollapsedVert(mTriangulation.cells(2, cell));
            size_t d = getCollapsedVert(mTriangulation.cells(3, cell));

            if(a == b || a == c || a == d || b == c || b == d || c == d)
            {
                //the cell is collapsed
                continue;
            }

            matrixr_t invA = ones<double>(4, 4);
            invA(colon(0, 2), 0) = mTriangulation.vertices(colon(), a);
            invA(colon(0, 2), 1) = mTriangulation.vertices(colon(), b);
            invA(colon(0, 2), 2) = mTriangulation.vertices(colon(), c);
            invA(colon(0, 2), 3) = mTriangulation.vertices(colon(), d);

            if(inv(invA)) {
                std::cerr << "warning: degenerated tetrahedron." << std::endl;
                continue;
            }
            for(int i = 0; i < sampleInner.size(); i++)
            {
                const vec4_t bary = invA(colon(), colon(0, 2)) * mShell.mInnerShell(colon(), sampleInner[i]) + invA(colon(), 3);
                if(min(bary) >= 0)
                {
                    this->mtx.lock();
                    innerSample.insert(sampleInner[i]);
                    this->mtx.unlock();
                }
            }
            for(int i = 0; i < sampleOuter.size(); i++)
            {
                const vec4_t bary = invA(colon(), colon(0, 2)) * mShell.mOuterShell(colon(), sampleOuter[i]) + invA(colon(), 3);
                if(min(bary) >= 0)
                {
                    this->mtx.lock();
                    outerSample.insert(sampleOuter[i]);
                    this->mtx.unlock();
                }
            }
        }
    }

    bool EdgeCollapse::findCollapsePos(size_t vert, size_t vertCollapseTo, vec3_t &position, double& out_error)
    {
        matrixs_t faces;
        std::set<size_t> innerSample;
        std::set<size_t> outerSample;
        bool found = false;
        buildOneRingArea(vert, vertCollapseTo, faces, innerSample, outerSample);

        KernelRegion kernel(mTriangulation.vertices, faces, mTriangulation.vertices(colon(), vert), mShell, innerSample, outerSample,
                                 mTriangulation, mTriangulation.vertType[vert]);
        out_error = std::numeric_limits<double>::max();

        if(mType == BOUNDARY)
        {
            if(mTriangulation.vertType[vert] == POINT_INNER)
            {
                for(size_t sample : innerSample)
                {
                    const vec3_t samplePoint = mShell.mInnerShell(colon(), sample);
                    if(kernel.contains(samplePoint))
                    {
                        if(checkNormal(faces, samplePoint, POINT_INNER) == false)
                            continue;

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
                    const vec3_t samplePoint = mShell.mOuterShell(colon(), sample);
                    if(kernel.contains(samplePoint))
                    {
                        if(checkNormal(faces, samplePoint, POINT_OUTER) == false)
                            continue;

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
            double xmax = std::numeric_limits<double>::lowest();
            double ymin = std::numeric_limits<double>::max();
            double ymax = std::numeric_limits<double>::lowest();
            double zmin = std::numeric_limits<double>::max();
            double zmax = std::numeric_limits<double>::lowest();

            for(int i = 0; i < faces.size(2); i++)
            {
                const matrixr_t& a = mTriangulation.vertices(colon(), faces(0, i));
                const matrixr_t& b = mTriangulation.vertices(colon(), faces(1, i));
                const matrixr_t& c = mTriangulation.vertices(colon(), faces(2, i));
                xmax = std::max(std::max(a[0], b[0]), c[0]);
                xmin = std::min(std::min(a[0], b[0]), c[0]);
                ymax = std::max(std::max(a[1], b[1]), c[1]);
                ymin = std::min(std::min(a[1], b[1]), c[1]);
                zmax = std::max(std::max(a[2], b[2]), c[2]);
                zmin = std::min(std::min(a[2], b[2]), c[2]);
            }

            SamplingQuadTree tree(kernel, xmax, xmin, ymax, ymin, zmax, zmin, mSampleRadius);
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

    bool EdgeCollapse::checkNormal(const matrixs_t &oneRingFaces, const vec3_t &point, PointType pointType)
    {
        for(int i = 0; i < oneRingFaces.size(2); i++)
        {
            matrixr_t tetra(3, 4);
            tetra(colon(), colon(0, 2)) = mTriangulation.vertices(colon(), oneRingFaces(colon(), i));
            tetra(colon(), 3) = point;

            PointType type_v1 = mTriangulation.vertType[oneRingFaces(0, i)];
            PointType type_v2 = mTriangulation.vertType[oneRingFaces(1, i)];
            PointType type_v3 = mTriangulation.vertType[oneRingFaces(2, i)];
            PointType type_v4 = pointType;

            if(mTriangulation.getFValue(type_v1) == mTriangulation.getFValue(type_v2)
                    && mTriangulation.getFValue(type_v1) == mTriangulation.getFValue(type_v3)
                    && mTriangulation.getFValue(type_v1) == mTriangulation.getFValue(type_v4))
                continue;

            if(NormalChecker::check(tetra, type_v1, type_v2, type_v3, type_v4, mShell) == false)
            {
                return false;
            }
        }
        return true;
    }
}
