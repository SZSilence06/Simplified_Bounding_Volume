#include "KdTreeWrap.h"
#include <wkylib/geometry.h>
#include <limits>
#include <algorithm>   //for std::max()

using namespace zjucad::matrix;

namespace SBV
{
    void KdTreeWrap::build(const std::vector<Point> &points)
    {
        for(size_t i = 0; i < points.size(); i++)
        {
            KdTreeNode node;
            node.point = points[i];
            node.index = i;
            mTree.insert(node);
        }
    }

    size_t KdTreeWrap::getNearestPoint(const Point &point) const
    {
        KdTreeNode node;
        node.point = point;
        auto result = mTree.find_nearest(node);
        return result.first->index;
    }

#ifdef VER_2D
    void KdTreeWrap::getPointsInRange(double xmin, double xmax, double ymin, double ymax, std::vector<size_t> &points) const
    {
        KdTreeNode node;
        node.point[0] = (xmax + xmin) / 2;
        node.point[1] = (ymax + ymin) / 2;
        double range = std::max(node[0] - xmin, node[1] - ymin);

        std::vector<KdTreeNode> find_results(mTree.size());
        auto end = mTree.find_within_range(node, range, find_results.begin());

        //organize output
        points.resize(end - find_results.begin());
        for(size_t i = 0; i < points.size(); i++)
        {
            points[i] = find_results[i].index;
        }
    }
#else
    void KdTreeWrap::getPointsInRange(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                                      std::vector<size_t>& points) const
    {
        KdTreeNode node;
        node.point[0] = (xmax + xmin) / 2;
        node.point[1] = (ymax + ymin) / 2;
        node.point[2] = (ymax + ymin) / 2;
        double range = std::max(std::max(node[0] - xmin, node[1] - ymin), node[2] - zmin);

        std::vector<KdTreeNode> find_results(mTree.size());
        auto end = mTree.find_within_range(node, range, find_results.begin());

        //organize output
        points.resize(end - find_results.begin());
        for(size_t i = 0; i < points.size(); i++)
        {
            points[i] = find_results[i].index;
        }
    }
#endif
}
