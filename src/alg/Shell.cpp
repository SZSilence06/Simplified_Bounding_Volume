#include "Shell.h"

namespace SBV
{
    void Shell::buildKdTree()
    {
        mInnerTree.build(mInnerShell);
        mOuterTree.build(mOuterShell);
    }
}
