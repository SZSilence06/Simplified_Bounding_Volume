#ifndef WKY_BARY_COMPUTER_H
#define WKY_BARY_COMPUTER_H

#include "Common.h"
#include <set>

namespace SBV
{
    class Shell;

    class BaryComputer
    {
    public:
      BaryComputer(const matrixr_t& triangle);
      template <typename M>
      inline void operator()(const M& pt2d, vec3_t &bary) const {
        using namespace zjucad::matrix;
        bary = invA(colon(), colon(0, 1))*pt2d + invA(colon(), 2);
      }

    private:
        void buildInvA(const matrixr_t& triangle);
    private:
        matrixr_t invA;
    };
}

#endif
