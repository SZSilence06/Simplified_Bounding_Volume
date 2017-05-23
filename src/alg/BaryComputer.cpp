#include "BaryComputer.h"
#include "Shell.h"
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <iostream>

using namespace zjucad::matrix;

namespace SBV
{
    BaryComputer::BaryComputer(const Shell& shell, const matrixr_t &triangle, const matrixs_t &sampleInner, const matrixs_t &sampleOuter)
        : mShell(shell)
    {
        buildInvA(triangle);
        this->sampleInnerPtr = &sampleInner;
        this->sampleOuterPtr = &sampleOuter;
    }

    BaryComputer::BaryComputer(const Shell& shell, const matrixr_t &triangle, const std::set<size_t> &sampleInner, const std::set<size_t> &sampleOuter)
        : mShell(shell)
    {
        buildInvA(triangle);
        this->_sampleInner.resize(sampleInner.size());
        this->_sampleOuter.resize(sampleOuter.size());
        int i = 0;
        for(size_t sample : sampleInner)
        {
            this->_sampleInner[i] = sample;
            i++;
        }
        i = 0;
        for(size_t sample : sampleOuter)
        {
            this->_sampleOuter[i] = sample;
            i++;
        }
        this->sampleInnerPtr = &this->_sampleInner;
        this->sampleOuterPtr = &this->_sampleOuter;
    }

    void BaryComputer::computeBary(matrixr_t &barysInner, matrixr_t &barysOuter)
    {
        _computeBary(*sampleInnerPtr, barysInner, true);
        _computeBary(*sampleOuterPtr, barysOuter, false);
    }

    void BaryComputer::_computeBary(const matrixs_t &samples, matrixr_t &barys, bool isInner)
    {
        if(isInner)
        {
            barys = invA(colon(), colon(0, 1)) * mShell.mInnerShell(colon(), samples) +
                    invA(colon(), 2) * ones<double>(1, samples.size());
        }
        else
        {
            barys = invA(colon(), colon(0, 1)) * mShell.mOuterShell(colon(), samples) +
                    invA(colon(), 2) * ones<double>(1, samples.size());
        }
    }

    void BaryComputer::buildInvA(const matrixr_t &triangle)
    {
        this->invA = ones<double>(3, 3);
        this->invA(colon(0, 1), colon()) = triangle;
        if(inv(invA)) {
            std::cerr << "warning: degenerated triangle." << std::endl;
        }
    }
}
