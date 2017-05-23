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
        BaryComputer(const Shell& shell, const matrixr_t& triangle, const matrixs_t& sampleInner, const matrixs_t& sampleOuter);
        BaryComputer(const Shell& shell, const matrixr_t &triangle, const std::set<size_t> &sampleInner, const std::set<size_t> &sampleOuter);

        void computeBary(matrixr_t& barysInner, matrixr_t& barysOuter);

    private:
        void _computeBary(const matrixs_t& samples, matrixr_t& barys, bool isInner);
        void buildInvA(const matrixr_t& triangle);

    private:
        const Shell& mShell;
        matrixr_t invA;
        double invA_data[9];
        const matrixs_t* sampleInnerPtr = nullptr;
        const matrixs_t* sampleOuterPtr = nullptr;

        //internal stored samples
        matrixs_t _sampleInner;
        matrixs_t _sampleOuter;
    };
}

#endif
