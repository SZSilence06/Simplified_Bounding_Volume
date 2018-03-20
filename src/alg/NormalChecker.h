#ifndef WKY_NORMAL_CHECKER_H
#define WKY_NORMAL_CHECKER_H

#include "Common.h"
#include "TriangulatedShell.h"

namespace SBV
{
    class Shell;

    /**
     * @brief The NormalChecker class
     *
     * This class is used to check whether a cell satisfies the normal condition.
     */
    class NormalChecker
    {
    public:
        /**
         * @brief check whether a cell satisfies the normal condition.
         * @param cell : 3x4 matrix representing the cell.
         * @param type_v1 : point type of the 1st vertex of the cell.
         * @param type_v2 : point type of the 2nd vertex of the cell.
         * @param type_v3 : point type of the 3rd vertex of the cell.
         * @param type_v4 : point type of the 4th vertex of the cell.
         * @param shell : input samples.
         * @return true if satisfies, and false otherwise.
         */
        static bool check(const matrixr_t& cell, PointType type_v1, PointType type_v2, PointType type_v3, PointType type_v4,
                                    const Shell& shell);
    };
}

#endif
