#include "KdTreeWrap.h"
#include <wkylib/geometry.h>
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

    void KdTreeWrap::getPointsInPolygon(const matrixr_t &polygon, matrixs_t &points) const
    {
        double xmax = std::numeric_limits<double>::min();
        double xmin = std::numeric_limits<double>::max();
        double ymax = std::numeric_limits<double>::min();
        double ymin = std::numeric_limits<double>::max();

        //build bounding box
        for(int i = 0; i < polygon.size(2); i++)
        {
            const matrixr_t& vert = polygon(colon(), i);
            if(vert[0] > xmax)
            {
                xmax = vert[0];
            }
            if(vert[0] < xmin)
            {
                xmin = vert[0];
            }
            if(vert[1] > ymax)
            {
                ymax = vert[1];
            }
            if(vert[1] < ymin)
            {
                ymin = vert[1];
            }
        }

        KdTreeNode node;
        node.point = matrixr_t(2,1);
        node.point[0] = (xmax + xmin) / 2;
        node.point[1] = (ymax + ymin) / 2;
        double range = std::max(node[0] - xmin, node[1] - ymin);

        std::vector<KdTreeNode> find_results(mTree.size());
        auto end = mTree.find_within_range(node, range, find_results.begin());

        std::vector<size_t> indices;

        int i = 0;
        for(auto iter = find_results.begin(); iter != end; ++iter)
        {
            if(WKYLIB::is_inside_poly(polygon, iter->point))
            {
                indices.push_back(i);
            }
            i++;
        }

        //organize output
        points.resize(indices.size(), 1);
        for(int i = 0; i < indices.size(); i++)
        {
            points[i] = find_results[indices[i]].index;
        }
    }

    void KdTreeWrap::getPointsInRange(double xmin, double xmax, double ymin, double ymax, std::list<size_t> &points) const
    {
        KdTreeNode node;
        node.point = matrixr_t(2,1);
        node.point[0] = (xmax + xmin) / 2;
        node.point[1] = (ymax + ymin) / 2;
        double range = std::max(node[0] - xmin, node[1] - ymin);

        std::vector<KdTreeNode> find_results(mTree.size());
        auto end = mTree.find_within_range(node, range, find_results.begin());

        //organize output       
        points.clear();
        for(int i = 0; i < end - find_results.begin(); i++)
        {
            points.push_back(find_results[i].index);
        }
    }
}
