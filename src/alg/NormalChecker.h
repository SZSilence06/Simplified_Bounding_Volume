#ifndef WKY_NORMAL_CHECKER_H
#define WKY_NORMAL_CHECKER_H

#include "Common.h"
#include "TriangulatedShell.h"

namespace SBV
{
    class Shell;

    class NormalChecker
    {
    public:
        static bool check(const matrixr_t& cell, PointType type_v1, PointType type_v2, PointType type_v3, PointType type_v4,
                                    const Shell& shell);
    };
}

#endif
