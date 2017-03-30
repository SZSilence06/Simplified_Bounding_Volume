#include "Shell.h"

using namespace zjucad::matrix;

namespace SBV
{
    void Shell::buildKdTree()
    {
        mInnerTree.build(mInnerShell);
        mOuterTree.build(mOuterShell);
    }
}
