#include "TriangulatedShell.h"
#include <jtflib/mesh/io.h>

using namespace zjucad::matrix;

namespace SBV
{
    double TriangulatedShell::getFValue(PointType pointType)
    {
        switch(pointType)
        {
        case POINT_BOUNDING_BOX:
            return 1;
        case POINT_OUTER:
            return 1;
        case POINT_INNER:
            return -1;
        case POINT_ZERO:
            return 0;
        default:
            throw std::runtime_error("not on refined shell");
        }
    }

    double TriangulatedShell::getFValue(size_t vert) const
    {
        return getFValue(vertType[vert]);
    }

    double TriangulatedShell::getSign(size_t vert) const
    {
        double value = getFValue(vert);
        if(value > 0)
        {
            return 1;
        }
        else if(value == 0)
        {
            return 0;
        }
        return -1;
    }

    void TriangulatedShell::buildZeroSet()
    {
        //find out the number of faces in the zero set
        if(hasZeroSet)
        {
           buildZeroSetExisting();
           return;
        }
        mZeroSet.vertPairs.clear();
        mZeroSet.faceTetras.clear();

        //build zero face connection
        std::vector<matrixs_t> triangleVector;
        int test = 0;
        for(int i = 0; i < cells.size(2); i++)
        {
            size_t v0 = cells(0, i);
            size_t v1 = cells(1, i);
            size_t v2 = cells(2, i);
            size_t v3 = cells(3, i);

            if(getSign(v0) == getSign(v1) && getSign(v0) == getSign(v2) && getSign(v0) == getSign(v3))
            {
                //F value sign of the vertices are same, so no zero-set in this triangle
                continue;
            }

            std::vector<size_t> innerVerts, outerVerts;
            for(int j = 0; j < 4; j++)
            {
                getSign(cells(j, i)) < 0 ? innerVerts.push_back(cells(j, i)) : outerVerts.push_back(cells(j, i));
            }

            if(innerVerts.size() == 2 && outerVerts.size() == 2)
            {
                matrixs_t tri1(3, 1), tri2(3, 1);
                const vec3_t a = (vertices(colon(), innerVerts[0]) + vertices(colon(), outerVerts[0])) / 2;
                const vec3_t b = (vertices(colon(), innerVerts[0]) + vertices(colon(), outerVerts[1])) / 2;
                const vec3_t c = (vertices(colon(), innerVerts[1]) + vertices(colon(), outerVerts[1])) / 2;
                const vec3_t d = (vertices(colon(), innerVerts[1]) + vertices(colon(), outerVerts[0])) / 2;

                if(norm(a - c) < norm(b - d))
                {
                    tri1[0] = getZeroPointIndex(innerVerts[0], outerVerts[0]);
                    tri1[1] = getZeroPointIndex(innerVerts[0], outerVerts[1]);
                    tri1[2] = getZeroPointIndex(innerVerts[1], outerVerts[1]);
                    tri2[0] = getZeroPointIndex(innerVerts[0], outerVerts[0]);
                    tri2[1] = getZeroPointIndex(innerVerts[1], outerVerts[1]);
                    tri2[2] = getZeroPointIndex(innerVerts[1], outerVerts[0]);
                }
                else
                {
                    tri1[0] = getZeroPointIndex(innerVerts[0], outerVerts[0]);
                    tri1[1] = getZeroPointIndex(innerVerts[0], outerVerts[1]);
                    tri1[2] = getZeroPointIndex(innerVerts[1], outerVerts[0]);
                    tri2[0] = getZeroPointIndex(innerVerts[0], outerVerts[1]);
                    tri2[1] = getZeroPointIndex(innerVerts[1], outerVerts[0]);
                    tri2[2] = getZeroPointIndex(innerVerts[1], outerVerts[1]);
                }
                triangleVector.push_back(tri1);
                triangleVector.push_back(tri2);
                if(test == 529 || test == 530)
                {
                    std::cout << "cell id " << i << std::endl;
                }
                test += 2;
            }
            else
            {
                matrixs_t tri(3, 1);
                if(innerVerts.size() == 3)
                {
                    tri[0] = getZeroPointIndex(innerVerts[0], outerVerts[0]);
                    tri[1] = getZeroPointIndex(innerVerts[1], outerVerts[0]);
                    tri[2] = getZeroPointIndex(innerVerts[2], outerVerts[0]);
                }
                else
                {
                    tri[0] = getZeroPointIndex(innerVerts[0], outerVerts[0]);
                    tri[1] = getZeroPointIndex(innerVerts[0], outerVerts[1]);
                    tri[2] = getZeroPointIndex(innerVerts[0], outerVerts[2]);
                }
                triangleVector.push_back(tri);
                if(test == 529)
                {
                    std::cout << "cell id " << i << std::endl;
                }
                test++;
            }
        }

        mZeroSet.triangles.resize(3, triangleVector.size());
        for(int i = 0; i < triangleVector.size(); i++)
        {
            mZeroSet.triangles(colon(), i) = triangleVector[i];
        }

        //build zero point positions
        mZeroSet.vertices.resize(3, mZeroSet.vertPairs.size());
        for(int i = 0; i< mZeroSet.vertPairs.size(); i++)
        {
            auto& vertPair = mZeroSet.vertPairs[i];
            mZeroSet.vertices(colon(), i) = (vertices(colon(), vertPair.first) + vertices(colon(), vertPair.second)) / 2;
        }
    }

    void TriangulatedShell::buildZeroSetExisting()
    {
        mZeroSet.faceTetras.clear();

        std::vector<matrixr_t> zeroVerts;
        for(int i = 0; i < vertices.size(2); i++)
        {
            if(vertType[i] == POINT_ZERO)
            {
                zeroVerts.push_back(vertices(colon(), i));
            }
        }

        //organize zero vertices
        matrixr_t newVerts(3, zeroVerts.size());
        for(int i = 0; i < zeroVerts.size(); i++)
        {
            newVerts(colon(), i) = zeroVerts[i];
        }
        mZeroSet.vertices = newVerts;

        //find out zero faces count
        std::set<ZeroFace > zeroFaces;
        for(int i = 0; i < cells.size(2); i++)
        {
            std::vector<size_t> zeroVert;
            for(int j = 0; j < 4; j++)
            {
                if(vertType[cells(j, i)] == POINT_ZERO)
                {
                    zeroVert.push_back(cells(j, i));
                }
            }
            if(zeroVert.size() == 3)
            {
                tryAddZeroFace(i, zeroVert[0], zeroVert[1], zeroVert[2], zeroFaces);
            }
        }

        //organize zero faces
        int i = 0;
        mZeroSet.triangles.resize(3, zeroFaces.size());
        for(const ZeroFace& face : zeroFaces)
        {
            int j = 0;
            for(const size_t& vert : face.verts)
            {
                mZeroSet.triangles(j, i) = vert;
                j++;
            }
            mZeroSet.faceTetras.push_back(face.tetra);
            i++;
        }
    }

    void TriangulatedShell::tryAddZeroFace(size_t currentTetra, size_t zeroVert1, size_t zeroVert2, size_t zeroVert3,
                                           std::set<ZeroFace>& zeroFaces)
    {
        int zeroIndexDelta = vertices.size(2) - mZeroSet.vertices.size(2);
        ZeroFace face;
        face.verts.insert(zeroVert1 - zeroIndexDelta);
        face.verts.insert(zeroVert2 - zeroIndexDelta);
        face.verts.insert(zeroVert3 - zeroIndexDelta);
        face.tetra = currentTetra;
        zeroFaces.insert(face);
    }

    size_t TriangulatedShell::getZeroPointIndex(size_t firstVertex, size_t secondVertex)
    {
        for(int i = 0; i < mZeroSet.vertPairs.size(); i++)
        {
            auto& vertPair = mZeroSet.vertPairs[i];
            if((vertPair.first == firstVertex && vertPair.second == secondVertex)
                    || (vertPair.first == secondVertex && vertPair.second == firstVertex))
            {
                return i;
            }
        }

        mZeroSet.vertPairs.push_back(std::make_pair(firstVertex, secondVertex));
        return mZeroSet.vertPairs.size() - 1;
    }

    void TriangulatedShell::mutualTessellate()
    {
        hasZeroSet = true;
        std::vector<matrixs_t> newCells;
        for(int i = 0; i < cells.size(2); i++)
        {
            const matrixs_t& cell = cells(colon(), i);
            if(getSign(cell[0]) == getSign(cell[1]) && getSign(cell[0]) == getSign(cell[2]) && getSign(cell[0]) == getSign(cell[3]))
            {
                //no zero set, so this face won't be splited
                newCells.push_back(cell);
                continue;
            }

            std::vector<size_t> innerVerts, outerVerts;
            for(int j = 0; j < 4; j++)
            {
                getSign(cells(j, i)) < 0 ? innerVerts.push_back(cells(j, i)) : outerVerts.push_back(cells(j, i));
            }
            if(innerVerts.size() == 2)
            {
                size_t a = innerVerts[0], b = innerVerts[1];
                size_t c = outerVerts[0], d = outerVerts[1];
                size_t e = getZeroPointIndex(a, c) + vertices.size(2);
                size_t f = getZeroPointIndex(a, d) + vertices.size(2);
                size_t g = getZeroPointIndex(b, c) + vertices.size(2);
                size_t h = getZeroPointIndex(b, d) + vertices.size(2);

                matrixs_t cell1(4, 1), cell2(4, 1), cell3(4, 1), cell4(4, 1), cell5(4, 1), cell6(4, 1);
                cell1[0] = a;
                cell1[1] = b;
                cell1[2] = e;
                cell1[3] = f;

                cell2[0] = b;
                cell2[1] = e;
                cell2[2] = f;
                cell2[3] = h;

                cell3[0] = d;
                cell3[1] = e;
                cell3[2] = f;
                cell3[3] = h;

                cell4[0] = c;
                cell4[1] = d;
                cell4[2] = e;
                cell4[3] = g;

                cell5[0] = b;
                cell5[1] = e;
                cell5[2] = g;
                cell5[3] = h;

                cell6[0] = d;
                cell6[1] = e;
                cell6[2] = g;
                cell6[3] = h;

                newCells.push_back(cell1);
                newCells.push_back(cell2);
                newCells.push_back(cell3);
                newCells.push_back(cell4);
                newCells.push_back(cell5);
                newCells.push_back(cell6);
            }
            else
            {
                size_t a, b, c, d;
                if(innerVerts.size() == 1)
                {
                    a = innerVerts[0];
                    b = outerVerts[0];
                    c = outerVerts[1];
                    d = outerVerts[2];
                }
                else
                {
                    a = outerVerts[0];
                    b = innerVerts[0];
                    c = innerVerts[1];
                    d = innerVerts[2];
                }
                size_t e = getZeroPointIndex(a, b) + vertices.size(2);
                size_t f = getZeroPointIndex(a, c) + vertices.size(2);
                size_t g = getZeroPointIndex(a, d) + vertices.size(2);

                matrixs_t cell1(4, 1), cell2(4, 1), cell3(4, 1), cell4(4, 1);
                cell1[0] = a;
                cell1[1] = e;
                cell1[2] = f;
                cell1[3] = g;

                cell2[0] = d;
                cell2[1] = e;
                cell2[2] = f;
                cell2[3] = g;

                cell3[0] = b;
                cell3[1] = c;
                cell3[2] = d;
                cell3[3] = e;

                cell4[0] = c;
                cell4[1] = d;
                cell4[2] = e;
                cell4[3] = f;

                newCells.push_back(cell1);
                newCells.push_back(cell2);
                newCells.push_back(cell3);
                newCells.push_back(cell4);
            }
        }

        //organize output
        matrixr_t newVerts(3, vertices.size(2) + mZeroSet.vertices.size(2));
        newVerts(colon(), colon(0, vertices.size(2) - 1)) = vertices;
        newVerts(colon(), colon(vertices.size(2), newVerts.size(2) - 1)) = mZeroSet.vertices;
        vertices = newVerts;

        cells.resize(4, newCells.size());
        for(int i = 0; i < newCells.size(); i++)
        {
            cells(colon(), i) = newCells[i];
        }

        //update vertType
        for(int i = 0; i < mZeroSet.vertices.size(2); i++)
        {
            vertType.push_back(POINT_ZERO);
        }
    }

    void TriangulatedShell::outputBoundaryFaces(const std::string &path)
    {
        std::vector<matrixs_t> innerBoundary;
        std::vector<matrixs_t> outerBoundary;

        for(int i = 0; i < cells.size(2); i++)
        {
            std::vector<size_t> innerVerts, outerVerts;
            for(int j = 0; j < 4; j++)
            {
                getSign(cells(j, i)) < 0 ? innerVerts.push_back(cells(j, i)) : outerVerts.push_back(cells(j, i));
            }

            if(innerVerts.size() == 3)
            {
                matrixs_t face(3, 1);
                face[0] = innerVerts[0];
                face[1] = innerVerts[1];
                face[2] = innerVerts[2];
                innerBoundary.push_back(face);
            }
            else if(outerVerts.size() == 3)
            {
                matrixs_t face(3, 1);
                face[0] = outerVerts[0];
                face[1] = outerVerts[1];
                face[2] = outerVerts[2];
                outerBoundary.push_back(face);
            }
        }

        matrixs_t inner(3, innerBoundary.size());
        matrixs_t outer(3, outerBoundary.size());
        for(int i = 0; i < innerBoundary.size(); i++)
        {
            inner(colon(), i) = innerBoundary[i];
        }
        for(int i = 0; i < outerBoundary.size(); i++)
        {
            outer(colon(), i) = outerBoundary[i];
        }

        jtf::mesh::save_obj((path + "_innerBoundary.obj").c_str(), inner, this->vertices);
        jtf::mesh::save_obj((path + "_outerBoundary.obj").c_str(), outer, this->vertices);
    }
}
