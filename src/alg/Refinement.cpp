#include "Refinement.h"
#include <limits>

namespace SBV
{
    Refinement::Refinement(const matrixr_t &innerShell, const matrixr_t &outerShell)
        : mInnerShell(innerShell),
          mOuterShell(outerShell)
    {
        init();
    }

    void Refinement::init()
    {
        computeBoundingBox();
        initErrors();
    }

    void Refinement::computeBoundingBox()
    {
        double xMax = std::numeric_limits<double>::min();
        double xMin = std::numeric_limits<double>::max();
        double yMax = std::numeric_limits<double>::min();
        double yMin = std::numeric_limits<double>::max();
        double zMax = std::numeric_limits<double>::min();
        double zMin = std::numeric_limits<double>::max();

        for(int i = 0; i < mOuterShell.size(2); i++)
        {
            if(xMax < mOuterShell(0, i))
            {
                xMax = mOuterShell(0, i);
            }
            if(xMin > mOuterShell(0, i))
            {
                xMin = mOuterShell(0, i);
            }
            if(yMax < mOuterShell(1, i))
            {
                yMax = mOuterShell(1, i);
            }
            if(yMin > mOuterShell(1, i))
            {
                yMin = mOuterShell(1, i);
            }
            if(zMax < mOuterShell(2, i))
            {
                zMax = mOuterShell(2, i);
            }
            if(zMin > mOuterShell(2, i))
            {
                zMin = mOuterShell(2, i);
            }
        }

        mDelaunay.insert(Point(xMax, yMax, zMax));
        mDelaunay.insert(Point(xMax, yMax, zMin));
        mDelaunay.insert(Point(xMax, yMin, zMax));
        mDelaunay.insert(Point(xMax, yMin, zMin));
        mDelaunay.insert(Point(xMin, yMax, zMax));
        mDelaunay.insert(Point(xMin, yMax, zMin));
        mDelaunay.insert(Point(xMin, yMin, zMax));
        mDelaunay.insert(Point(xMin, yMin, zMin));
    }

    void Refinement::initErrors()
    {
        mInnerError.reserve(mInnerShell.size(2));
        for(int i = 0; i < mInnerShell.size(2); i++)
        {
            innerError.push_back(0);
        }

        mOuterError.reserve(mOuterShell.size(2));
        for(int i = 0; i < mOuterShell.size(2); i++)
        {
            mOuterError.push_back(0);
        }
    }

    bool Refinement::refine(std::vector<size_t> &output_refinement)
    {
        std::vector<double> innerError;
        std::vector<double> outerError;

        innerError.reserve(mInnerShell.size(2));
        for(int i = 0; i < mInnerShell.size(2); i++)
        {
            innerError.push_back(0);
        }

        outerError.reserve(mOuterShell.size(2));
        for(int i = 0; i < mOuterShell.size(2); i++)
        {
            outerError.push_back(0);
        }

        while(!isFinished())
        {

        }

        return true;
    }

    bool Refinement::isFinished()
    {
        return false;
    }
}
