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
#ifdef VER_2D
        struct ZeroSet
        {
            std::vector<Point> vertices;
            std::vector<Eigen::Vector2i> lines;
            std::vector<std::pair<size_t, size_t> >  vertPairs;   //recording the corresponding vert pairs of the zero set vertices
            std::vector<size_t> lineFaces;                        //recording the corresponding faces of the zero set lines
        };
#else
        struct ZeroSet
        {
            std::vector<Point> vertices;
            std::vector<Eigen::Vector3i> triangles;
            std::vector<std::pair<size_t, size_t> >  vertPairs;   //recording the corresponding vert pairs of the zero set vertices
            std::vector<size_t> faceTetras;                        //recording the corresponding faces of the zero set lines
        };
#endif

    public:
        std::vector<Point> vertices;
#ifdef VER_2D
        std::vector<Eigen::Vector3i> triangles;
#else
        std::vector<Eigen::Vector4i> cells;
#endif
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

#ifdef VER_2D
        void tryAddZeroFace(size_t currentTetra, size_t zeroVert1, size_t zeroVert2, std::set<ZeroFace>& zeroFaces);
#else
        void tryAddZeroFace(size_t currentTetra, size_t zeroVert1, size_t zeroVert2, size_t zeroVert3, std::set<ZeroFace>& zeroFaces);
#endif

    private:
        ZeroSet mZeroSet;
        bool hasZeroSet = false;
    };
}

#endif
