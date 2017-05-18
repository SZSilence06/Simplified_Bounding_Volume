#ifndef WKY_TRIANGULATED_SHELL_H
#define WKY_TRIANGULATED_SHELL_H

#include "Common.h"
#include "PointType.h"
#include <set>

namespace SBV
{
    class TriangulatedShell
    {
    public:
        struct ZeroSet
        {
            std::vector<Point> vertices;
            std::vector<Eigen::Vector2i> lines;
            std::vector<std::pair<size_t, size_t> >  vertPairs;   //recording the corresponding vert pairs of the zero set vertices
            std::vector<size_t> lineFaces;                        //recording the corresponding faces of the zero set lines
        };

    public:
        std::vector<Point> vertices;
        std::vector<Eigen::Vector3i> triangles;
        std::vector<PointType> vertType;

    public:
        double getFValue(PointType pointType) const;
        double getFValue(size_t vert) const;
        double getSign(size_t vert) const;

        inline const ZeroSet& getZeroSet() const { return mZeroSet; }

        void buildZeroSet();
        void mutualTessellate();

    private:
        struct ZeroFace{
            std::set<size_t> verts;
            size_t tetra;

            bool operator <(const ZeroFace& face) const { return this->verts < face.verts; }
        };

        size_t getZeroPointIndex(size_t firstVertex, size_t secondVertex);
        void buildZeroSetExisting();
        void tryAddZeroFace(size_t currentTetra, size_t zeroVert1, size_t zeroVert2, std::set<ZeroFace>& zeroFaces);

    private:
        ZeroSet mZeroSet;
        bool hasZeroSet = false;
    };
}

#endif
