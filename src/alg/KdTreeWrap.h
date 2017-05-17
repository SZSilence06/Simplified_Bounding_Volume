#ifndef WKY_KDTREE_WRAP_H
#define WKY_KDTREE_WRAP_H

#include "Common.h"
#include <cmath>
#include <kdtree++/kdtree.hpp>

namespace SBV
{
    struct KdTreeNode
    {
        typedef double value_type;

        Point point;
        size_t index;

        value_type operator[](size_t n) const
        {
            return point[n];
        }

        double distance( const KdTreeNode &node)
        {
#ifdef  VER_2D
            double x = point[0] - node.point[0];
            double y = point[1] - node.point[1];

    // this is not correct   return sqrt( x*x+y*y+z*z);

    // this is what kdtree checks with find_within_range()
    // the "manhattan distance" from the search point.
    // effectively, distance is the maximum distance in any one dimension.
            return std::max(fabs(x), fabs(y));
#else
            double x = point[0] - node.point[0];
            double y = point[1] - node.point[1];
            double z = point[2] - node.point[2];
            return std::max(std::max(fabs(x), fabs(y)), fabs(z));
#endif
        }
    };

#ifdef VER_2D
    using KdTreeType = KDTree::KDTree<2, KdTreeNode>;
#else
    using KdTreeType = KDTree::KDTree<3, KdTreeNode>;
#endif
    using Region = KdTreeType::_Region_;

    class KdTreeWrap
    {
    public:
        void build(const std::vector<Point>& points);

        size_t getNearestPoint(const Point& point) const;

#ifdef VER_2D
        void getPointsInRange(double xmin, double xmax, double ymin, double ymax, std::vector<size_t>& points) const;
#else
        void getPointsInRange(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                              std::vector<size_t>& points) const;
#endif

    private:
        KdTreeType mTree;
    };
}

#endif
