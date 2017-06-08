#include "TriangulatedShell.h"

using namespace zjucad::matrix;

namespace SBV
{
    double TriangulatedShell::getFValue(PointType pointType) const
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
                tri1[0] = getZeroPointIndex(innerVerts[0], outerVerts[0]);
                tri1[1] = getZeroPointIndex(innerVerts[0], outerVerts[1]);
                tri1[2] = getZeroPointIndex(innerVerts[1], outerVerts[1]);
                tri2[0] = getZeroPointIndex(innerVerts[0], outerVerts[0]);
                tri2[1] = getZeroPointIndex(innerVerts[1], outerVerts[1]);
                tri2[2] = getZeroPointIndex(innerVerts[1], outerVerts[0]);
                triangleVector.push_back(tri1);
                triangleVector.push_back(tri2);
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
        /*mZeroSet.faceTetras.clear();

        std::vector<matrixr_t> zeroVerts;
        for(int i = 0; i < vertices.size(2); i++)
        {
            if(vertType[i] == POINT_ZERO)
            {
                zeroVerts.push_back(vertices(colon(), i));
            }
        }

        //organize zero vertices
        matrixr_t newVerts(2, zeroVerts.size());
        for(int i = 0; i < zeroVerts.size(); i++)
        {
            newVerts(colon(), i) = zeroVerts[i];
        }
        mZeroSet.vertices = newVerts;

        //find out zero faces count
        std::set<ZeroFace > zeroFaces;
        for(int i = 0, j = 0; i < triangles.size(2); i++)
        {
            size_t v0 = triangles(0, i);
            size_t v1 = triangles(1, i);
            size_t v2 = triangles(2, i);

            if(getSign(v0) == getSign(v1) && getSign(v0) == 0)
            {
                tryAddZeroFace(i, v0, v1, zeroFaces);
            }
            else if(getSign(v0) == getSign(v2) && getSign(v0) == 0)
            {
                tryAddZeroFace(i, v0, v2, zeroFaces);
            }
            else if(getSign(v1) == getSign(v2) && getSign(v1) == 0)
            {
                tryAddZeroFace(i, v1, v2, zeroFaces);
            }
        }

        //organize zero faces
        int i = 0;
        mZeroSet.lines.resize(2, zeroFaces.size());
        for(const ZeroFace& face : zeroFaces)
        {
            int j = 0;
            for(const size_t& vert : face.verts)
            {
                mZeroSet.lines(j, i) = vert;
                j++;
            }
            mZeroSet.lineFaces.push_back(face.tetra);
            i++;
        }*/
    }

    void TriangulatedShell::tryAddZeroFace(size_t currentTetra, size_t zeroVert1, size_t zeroVert2,
                                           std::set<ZeroFace>& zeroFaces)
    {
        int zeroIndexDelta = vertices.size(2) - mZeroSet.vertices.size(2);
        ZeroFace face;
        face.verts.insert(zeroVert1 - zeroIndexDelta);
        face.verts.insert(zeroVert2 - zeroIndexDelta);
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
        /*hasZeroSet = true;
        std::vector<matrixs_t> newTriangles;
        for(int i = 0; i < triangles.size(2); i++)
        {
            const matrixs_t& face = triangles(colon(), i);
            if(getSign(face[0]) == getSign(face[1]) && getSign(face[0]) == getSign(face[2]))
            {
                //no zero set, so this face won't be splited
                newTriangles.push_back(face);
            }
        }

        for(int i = 0; i < mZeroSet.lineFaces.size(); i++)
        {
            const size_t& face = mZeroSet.lineFaces[i];
            const size_t& a = triangles(0, face);
            const size_t& b = triangles(1, face);
            const size_t& c = triangles(2, face);
            const size_t& lineVert11 = mZeroSet.vertPairs[mZeroSet.lines(0, i)].first;
            const size_t& lineVert12 = mZeroSet.vertPairs[mZeroSet.lines(0, i)].second;
            const size_t& lineVert21 = mZeroSet.vertPairs[mZeroSet.lines(1, i)].first;
            const size_t& lineVert22 = mZeroSet.vertPairs[mZeroSet.lines(1, i)].second;

            size_t commonVert;
            if(getSign(a) == getSign(b))
            {
                commonVert = c;
            }
            else if(getSign(a) == getSign(c))
            {
                commonVert = b;
            }
            else
            {
                commonVert = a;
            }

            //insert one new triangle
            matrixs_t newFace1(3, 1);
            newFace1[0] = vertices.size(2) + mZeroSet.lines(0, i);
            newFace1[1] = vertices.size(2) + mZeroSet.lines(1, i);
            newFace1[2] = commonVert;
            newTriangles.push_back(newFace1);

            std::vector<size_t> faceVerts;
            faceVerts.push_back(a);
            faceVerts.push_back(b);
            faceVerts.push_back(c);

            faceVerts.erase(std::remove(faceVerts.begin(), faceVerts.end(), commonVert), faceVerts.end());

            //insert another two new triangle
            matrixs_t newFace2(3, 1);
            matrixs_t newFace3(3, 1);
            newFace2[0] = vertices.size(2) + mZeroSet.lines(0, i);
            newFace2[1] = vertices.size(2) + mZeroSet.lines(1, i);
            newFace2[2] = faceVerts[0];
            newFace3[0] = faceVerts[0];
            newFace3[1] = faceVerts[1];
            if((lineVert11 == commonVert && lineVert12 == faceVerts[0]) ||
                    (lineVert12 == commonVert && lineVert11 == faceVerts[0]))
            {
                newFace3[2] = vertices.size(2) + mZeroSet.lines(1, i);
            }
            else
            {
                newFace3[2] = vertices.size(2) + mZeroSet.lines(0, i);
            }
            newTriangles.push_back(newFace2);
            newTriangles.push_back(newFace3);
        }

        //organize output
        matrixr_t newVerts(2, vertices.size(2) + mZeroSet.vertices.size(2));
        newVerts(colon(), colon(0, vertices.size(2) - 1)) = vertices;
        newVerts(colon(), colon(vertices.size(2), newVerts.size(2) - 1)) = mZeroSet.vertices;
        vertices = newVerts;

        triangles.resize(3, newTriangles.size());
        for(int i = 0; i < newTriangles.size(); i++)
        {
            triangles(colon(), i) = newTriangles[i];
        }

        //update vertType
        for(int i = 0; i < mZeroSet.vertices.size(2); i++)
        {
            vertType.push_back(POINT_ZERO);
        }*/
    }
}
