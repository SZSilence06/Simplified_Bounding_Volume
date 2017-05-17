#ifndef WKY_SHELL_H
#define WKY_SHELL_H

#include "Common.h"
#include "KdTreeWrap.h"

namespace SBV
{
    class Shell
    {
    public:
        std::vector<Point> mInnerShell;
        std::vector<Point> mOuterShell;

    public:
        void buildKdTree();

        inline const KdTreeWrap& getInnerTree() const { return mInnerTree; }
        inline const KdTreeWrap& getOuterTree() const { return mOuterTree; }

    private:
       KdTreeWrap mInnerTree;
       KdTreeWrap mOuterTree;
    };
}

#endif
