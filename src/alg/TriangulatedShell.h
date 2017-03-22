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

    struct TriangulatedShell
    {
        matrixr_t vertices;
        matrixs_t triangles;
        std::vector<PointType> vertType;

        double getFValue(PointType pointType) const
        {
            switch(pointType)
            {
            case POINT_BOUNDING_BOX:
                return 2;
            case POINT_OUTER:
                return 1;
            case POINT_INNER:
                return -1;
            default:
                throw std::runtime_error("not on refined shell");
            }
        }

        double getFValue(size_t vert) const
        {
            return getFValue(vertType[vert]);
        }

        double getSign(size_t vert) const
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
    };
}

#endif
