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

        vec3_t point;
        size_t index;
        value_type operator[](size_t n) const
        {
            return point[n];
        }

        double distance( const KdTreeNode &node)
        {
            double x = point[0] - node.point[0];
            double y = point[1] - node.point[1];
            double z = point[2] - node.point[2];

    // this is not correct   return sqrt( x*x+y*y+z*z);

    // this is what kdtree checks with find_within_range()
    // the "manhattan distance" from the search point.
    // effectively, distance is the maximum distance in any one dimension.
            return std::max(std::max(fabs(x), fabs(y)), fabs(z));
        }
    };

    using KdTreeType = KDTree::KDTree<3, KdTreeNode>;
    using Region = KdTreeType::_Region_;

    class KdTreeWrap
    {
    public:
        void build(const matrixr_t& points);

        size_t getNearestPoint(const matrixr_t& point) const;

        void getPointsInRange(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, matrixs_t& points) const;

    private:
        KdTreeType mTree;
    };
}

#endif
