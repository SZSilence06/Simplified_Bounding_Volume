#ifndef WKY_TRIANGULATED_SHELL_H
#define WKY_TRIANGULATED_SHELL_H

#include "Common.h"

namespace SBV
{
    enum PointType
    {
        POINT_BOUNDING_BOX,
        POINT_INNER,
        POINT_OUTER,
        POINT_UNKNOWN
    };

    class TriangulatedShell
    {
    public:
        matrixr_t vertices;
        matrixs_t triangles;
        std::vector<PointType> vertType;

    private:
        struct ZeroSet
        {
            matrixr_t vertices;
            matrixs_t lines;
            std::vector<std::pair<size_t, size_t> >  vertPairs;   //recording the corresponding vert pairs of the zero set vertices
            std::vector<size_t> lineFaces;                        //recording the corresponding faces of the zero set lines
        };

        ZeroSet mZeroSet;

        size_t getZeroPointIndex(size_t firstVertex, size_t secondVertex);

    public:
        double getFValue(PointType pointType) const;
        double getFValue(size_t vert) const;
        double getSign(size_t vert) const;

        inline const ZeroSet& getZeroSet() const { return mZeroSet; }

        void buildZeroSet();
        void mutualTessellate();
    };
}

#endif
