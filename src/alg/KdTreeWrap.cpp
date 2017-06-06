#include "KdTreeWrap.h"
#include <limits>
#include <algorithm>   //for std::max()

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

    size_t KdTreeWrap::getNearestPoint(const matrixr_t &point) const
    {
        KdTreeNode node;
        node.point = point;
        auto result = mTree.find_nearest(node);
        return result.first->index;
    }

    void KdTreeWrap::getPointsInRange(double xmin, double xmax, double ymin, double ymax, matrixs_t &points) const
    {
        KdTreeNode node;
        node.point = matrixr_t(2,1);
        node.point[0] = (xmax + xmin) / 2;
        node.point[1] = (ymax + ymin) / 2;
        double range = std::max(node[0] - xmin, node[1] - ymin);

        std::vector<KdTreeNode> find_results(mTree.size());
        auto end = mTree.find_within_range(node, range, find_results.begin());

        //organize output
        points.resize(end - find_results.begin(), 1);
        for(int i = 0; i < points.size(); i++)
        {
            points[i] = find_results[i].index;
        }
    }
}
