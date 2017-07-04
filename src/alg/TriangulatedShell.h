#ifndef WKY_TRIANGULATED_SHELL_H
#define WKY_TRIANGULATED_SHELL_H

#include "Common.h"
#include <set>

namespace SBV
{
    enum PointType
    {
        POINT_BOUNDING_BOX,
        POINT_INNER,
        POINT_OUTER,
        POINT_ZERO,
        POINT_UNKNOWN
    };

    class TriangulatedShell
    {
    public:
        struct ZeroSet
        {
            matrixr_t vertices;
            matrixs_t triangles;
            std::vector<std::pair<size_t, size_t> >  vertPairs;   //recording the corresponding vert pairs of the zero set vertices
            std::vector<size_t> faceTetras;                        //recording the corresponding tetras of the zero set triangles
        };

    public:
        matrixr_t vertices;
        matrixs_t cells;
        std::vector<PointType> vertType;

    public:
        double getFValue(PointType pointType) const;
        double getFValue(size_t vert) const;
        double getSign(size_t vert) const;

        inline const ZeroSet& getZeroSet() const { return mZeroSet; }

        void buildZeroSet();
        void mutualTessellate();
        void outputBoundaryFaces(const std::string& path);

    private:
        struct ZeroFace{
            std::set<size_t> verts;
            size_t tetra;

            bool operator <(const ZeroFace& face) const { return this->verts < face.verts; }
        };

        size_t getZeroPointIndex(size_t firstVertex, size_t secondVertex);
        void buildZeroSetExisting();
        void tryAddZeroFace(size_t currentTetra, size_t zeroVert1, size_t zeroVert2, size_t zeroVert3, std::set<ZeroFace>& zeroFaces);

    private:
        ZeroSet mZeroSet;
        bool hasZeroSet = false;
    };
}

#endif
