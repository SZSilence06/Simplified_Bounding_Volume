#include "KdTreeWrap.h"

using namespace zjucad::matrix;

namespace SBV
{
    void KdTreeWrap::build(const matrixr_t &points)
    {
        for(int i = 0; i < points.size(2); i++)
        {
            KdTreeNode node;
            node.point = points(colon(), i);
            node.index = i;
            mTree.insert(node);
        }
    }

    void KdTreeWrap::getNearestPoint(const matrixr_t &point, matrixr_t &nearest) const
    {
        KdTreeNode node;
        node.point = point;
        auto result = mTree.find_nearest(node);
        nearest = result.first->point;
    }
}
