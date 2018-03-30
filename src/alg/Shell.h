#ifndef WKY_SHELL_H
#define WKY_SHELL_H

#include "Common.h"
#include "KdTreeWrap.h"

namespace SBV
{
     /**
     * @brief The Shell class describing the input samples.
     *
     * The samples are organized as two point matrices. one represents the inner shell samples and another represents the outer shell samples.
     *
     * After assigning the two matrices, usually you need to call buildKdTree() to build the kd-tree acceleration structure.
     * The shell has two kd-trees, one for inner shell and another for outer shell.
     */
    class Shell
    {
    public:
        matrixr_t mInnerShell;
        matrixr_t mOuterShell;

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
