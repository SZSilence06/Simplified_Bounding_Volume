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
      BaryComputer(const matrixr_t& tetra);
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
