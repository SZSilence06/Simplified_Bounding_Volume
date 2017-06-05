#include "TriangulatedShell.h"

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

#ifdef VER_2D
    void TriangulatedShell::buildZeroSet()
    {
        //find out the number of faces in the zero set
        if(hasZeroSet)
        {
           buildZeroSetExisting();
           return;
        }
        mZeroSet.vertPairs.clear();
        mZeroSet.lineFaces.clear();

        //build zero face connection
        mZeroSet.lines.clear();
        for(size_t i = 0; i < triangles.size(); i++)
        {
            size_t v0 = triangles[i][0];
            size_t v1 = triangles[i][1];
            size_t v2 = triangles[i][2];

            if(getSign(v0) == getSign(v1) && getSign(v0) == getSign(v2))
            {
                //F value sign of the vertices are same, so no zero-set in this triangle
                continue;
            }

            Eigen::Vector2i line;
            if(getSign(v0) == getSign(v1))
            {
                line[0] = getZeroPointIndex(v0, v2);
                line[1] = getZeroPointIndex(v1, v2);
            }
            else if(getSign(v0) == getSign(v2))
            {
                line[0] = getZeroPointIndex(v0, v1);
                line[1] = getZeroPointIndex(v1, v2);
            }
            else if(getSign(v1) == getSign(v2))
            {
                line[0] = getZeroPointIndex(v0, v1);
                line[1] = getZeroPointIndex(v0, v2);
            }
            else
            {
                throw std::runtime_error("cannot find valid zero point when building zero sets.");
            }
            mZeroSet.lineFaces.push_back(i);
            mZeroSet.lines.push_back(line);
        }

        //build zero point positions
        mZeroSet.vertices.resize(mZeroSet.vertPairs.size());
        for(size_t i = 0; i< mZeroSet.vertPairs.size(); i++)
        {
            auto& vertPair = mZeroSet.vertPairs[i];
            mZeroSet.vertices[i] = (vertices[vertPair.first] + vertices[vertPair.second]) / 2;
        }
    }
#else
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
        mZeroSet.triangles.clear();
        for(size_t i = 0; i < cells.size(); i++)
        {
            size_t v0 = cells[i][0];
            size_t v1 = cells[i][1];
            size_t v2 = cells[i][2];
            size_t v3 = cells[i][3];

            if(getSign(v0) == getSign(v1) && getSign(v0) == getSign(v2) && getSign(v0) == getSign(v3))
            {
                //F value sign of the vertices are same, so no zero-set in this triangle
                continue;
            }

            std::vector<size_t> innerVerts, outerVerts;
            for(int j = 0; j < 4; j++)
            {
                getSign(cells[i][j]) < 0 ? innerVerts.push_back(cells[i][j]) : outerVerts.push_back(cells[i][j]);
            }
            if(innerVerts.size() == 2 && outerVerts.size() == 2)
            {
                Eigen::Vector3i face1, face2;
                face1[0] = getZeroPointIndex(innerVerts[0], outerVerts[0]);
                face1[1] = getZeroPointIndex(innerVerts[0], outerVerts[1]);
                face1[2] = getZeroPointIndex(innerVerts[1], outerVerts[1]);
                face2[0] = getZeroPointIndex(innerVerts[0], outerVerts[0]);
                face2[1] = getZeroPointIndex(innerVerts[1], outerVerts[1]);
                face2[2] = getZeroPointIndex(innerVerts[1], outerVerts[0]);
                mZeroSet.faceTetras.push_back(i);
                mZeroSet.faceTetras.push_back(i);
                mZeroSet.triangles.push_back(face1);
                mZeroSet.triangles.push_back(face2);
            }
            else
            {
                Eigen::Vector3i face;
                if(innerVerts.size() == 1)
                {
                    face[0] = getZeroPointIndex(innerVerts[0], outerVerts[0]);
                    face[1] = getZeroPointIndex(innerVerts[0], outerVerts[1]);
                    face[2] = getZeroPointIndex(innerVerts[0], outerVerts[2]);
                }
                else
                {
                    face[0] = getZeroPointIndex(outerVerts[0], innerVerts[0]);
                    face[1] = getZeroPointIndex(outerVerts[0], innerVerts[1]);
                    face[2] = getZeroPointIndex(outerVerts[0], innerVerts[2]);
                }
                mZeroSet.faceTetras.push_back(i);
                mZeroSet.triangles.push_back(face);
            }
        }

        //build zero point positions
        mZeroSet.vertices.resize(mZeroSet.vertPairs.size());
        for(size_t i = 0; i< mZeroSet.vertPairs.size(); i++)
        {
            auto& vertPair = mZeroSet.vertPairs[i];
            mZeroSet.vertices[i] = (vertices[vertPair.first] + vertices[vertPair.second]) / 2;
        }
    }
#endif

#ifdef VER_2D
    void TriangulatedShell::buildZeroSetExisting()
    {
        mZeroSet.lineFaces.clear();

        mZeroSet.vertices.clear();
        for(size_t i = 0; i < vertices.size(); i++)
        {
            if(vertType[i] == POINT_ZERO)
            {
                mZeroSet.vertices.push_back(vertices[i]);
            }
        }

        //find out zero faces count
        std::set<ZeroFace > zeroFaces;
        for(int i = 0, j = 0; i < triangles.size(); i++)
        {
            size_t v0 = triangles[i][0];
            size_t v1 = triangles[i][1];
            size_t v2 = triangles[i][2];

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
        mZeroSet.lines.resize(zeroFaces.size());
        for(const ZeroFace& face : zeroFaces)
        {
            int j = 0;
            for(const size_t& vert : face.verts)
            {
                mZeroSet.lines[i][j] = vert;
                j++;
            }
            mZeroSet.lineFaces.push_back(face.tetra);
            i++;
        }
    }
#else
    void TriangulatedShell::buildZeroSetExisting()
    {
        mZeroSet.faceTetras.clear();

        mZeroSet.vertices.clear();
        for(size_t i = 0; i < vertices.size(); i++)
        {
            if(vertType[i] == POINT_ZERO)
            {
                mZeroSet.vertices.push_back(vertices[i]);
            }
        }

        //find out zero faces count
        std::set<ZeroFace > zeroFaces;
        for(int i = 0; i < cells.size(); i++)
        {
            std::vector<size_t> zeroVerts, nonZeroVerts;
            for(int j = 0; j < 4; j++)
            {
                if(getSign(cells[i][j]) == 0)
                    zeroVerts.push_back(cells[i][j]);
            }
            if(zeroVerts.size() == 3)
            {
                tryAddZeroFace(i, zeroVerts[0], zeroVerts[1], zeroVerts[2], zeroFaces);
            }
        }
        //organize zero faces
        int i = 0;
        mZeroSet.triangles.resize(zeroFaces.size());
        for(const ZeroFace& face : zeroFaces)
        {
            int j = 0;
            for(const size_t& vert : face.verts)
            {
                mZeroSet.triangles[i][j] = vert;
                j++;
            }
            mZeroSet.faceTetras.push_back(face.tetra);
            i++;
        }
    }
#endif

#ifdef VER_2D
    void TriangulatedShell::tryAddZeroFace(size_t currentTetra, size_t zeroVert1, size_t zeroVert2,
                                        std::set<ZeroFace>& zeroFaces)
    {
        int zeroIndexDelta = vertices.size() - mZeroSet.vertices.size();
        ZeroFace face;
        face.verts.insert(zeroVert1 - zeroIndexDelta);
        face.verts.insert(zeroVert2 - zeroIndexDelta);
        face.tetra = currentTetra;
        zeroFaces.insert(face);
    }
#else
    void TriangulatedShell::tryAddZeroFace(size_t currentTetra, size_t zeroVert1, size_t zeroVert2, size_t zeroVert3,
                                        std::set<ZeroFace>& zeroFaces)
    {
        int zeroIndexDelta = vertices.size() - mZeroSet.vertices.size();
        ZeroFace face;
        face.verts.insert(zeroVert1 - zeroIndexDelta);
        face.verts.insert(zeroVert2 - zeroIndexDelta);
        face.verts.insert(zeroVert3 - zeroIndexDelta);
        face.tetra = currentTetra;
        zeroFaces.insert(face);
    }
#endif

    size_t TriangulatedShell::getZeroPointIndex(size_t firstVertex, size_t secondVertex)
    {
        for(size_t i = 0; i < mZeroSet.vertPairs.size(); i++)
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

    /*void TriangulatedShell::mutualTessellate()
    {
        this->hasZeroSet = true;
        std::vector<Eigen::Vector3i> newTriangles;
        for(size_t i = 0; i < triangles.size(); i++)
        {
            const Eigen::Vector3i& face = triangles[i];
            if(getSign(face[0]) == getSign(face[1]) && getSign(face[0]) == getSign(face[2]))
            {
                //no zero set, so this face won't be splited
                newTriangles.push_back(face);
            }
        }

        for(size_t i = 0; i < mZeroSet.lineFaces.size(); i++)
        {
            const size_t& face = mZeroSet.lineFaces[i];
            const size_t& a = triangles[face][0];
            const size_t& b = triangles[face][1];
            const size_t& c = triangles[face][2];
            const size_t& lineVert11 = mZeroSet.vertPairs[mZeroSet.lines[i][0]].first;
            const size_t& lineVert12 = mZeroSet.vertPairs[mZeroSet.lines[i][0]].second;

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
            Eigen::Vector3i newFace1;
            newFace1[0] = vertices.size() + mZeroSet.lines[i][0];
            newFace1[1] = vertices.size() + mZeroSet.lines[i][1];
            newFace1[2] = commonVert;
            newTriangles.push_back(newFace1);

            std::vector<size_t> faceVerts;
            faceVerts.push_back(a);
            faceVerts.push_back(b);
            faceVerts.push_back(c);
            faceVerts.erase(std::remove(faceVerts.begin(), faceVerts.end(), commonVert), faceVerts.end());

            //insert another two new triangle
            Eigen::Vector3i newFace2;
            Eigen::Vector3i newFace3;
            newFace2[0] = vertices.size() + mZeroSet.lines[i][0];
            newFace2[1] = vertices.size() + mZeroSet.lines[i][1];
            newFace2[2] = faceVerts[0];
            newFace3[0] = faceVerts[0];
            newFace3[1] = faceVerts[1];
            if((lineVert11 == commonVert && lineVert12 == faceVerts[0]) ||
                    (lineVert12 == commonVert && lineVert11 == faceVerts[0]))
            {
                newFace3[2] = verticesgit .size() + mZeroSet.lines[i][1];
            }
            else
            {
                newFace3[2] = vertices.size() + mZeroSet.lines[i][0];
            }
            newTriangles.push_back(newFace2);
            newTriangles.push_back(newFace3);
        }

        //organize output
        for(size_t i = 0; i < mZeroSet.vertices.size(); i++)
        {
            vertices.push_back(mZeroSet.vertices[i]);
        }
        triangles = std::move(newTriangles);

        //update vertType
        for(size_t i = 0; i < mZeroSet.vertices.size(); i++)
        {
            vertType.push_back(POINT_ZERO);
        }
    }*/
}
