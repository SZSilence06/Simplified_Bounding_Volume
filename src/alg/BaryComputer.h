#ifndef WKY_BARY_COMPUTER_H
#define WKY_BARY_COMPUTER_H

#include "Common.h"
#include <set>

namespace SBV
{
    class Shell;

    /**
     * @brief The BaryComputer class is used to compute the barycenter coordinate inside a tetrahedron.
     */
    class BaryComputer
    {
    public:
      /**
      * @brief BaryComputer constructor
      * @param tetra : the tetrahedron matrix of 3x4 size.
      */
      BaryComputer(const matrixr_t& tetra);

      /**
       * @brief operator to compute the barycenter coordinate
       * @param pt : the point to compute barycenter coordinate. can be matrix type or vec3_t type
       * @param bary : the output barycenter coordinate.
       */
      template <typename M>
      inline void operator()(const M& pt, vec4_t &bary) const {
        using namespace zjucad::matrix;
        bary = invA(colon(), colon(0, 2)) * pt + invA(colon(), 3);
      }

    private:
        void buildInvA(const matrixr_t& tetra);
    private:
        matrixr_t invA;
    };
}

#endif
