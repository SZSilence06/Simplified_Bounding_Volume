#ifndef WKY_SBV_COMMON_H
#define WKY_SBV_COMMON_H

#include <vector>
#include <eigen3/Eigen/Eigen>

//#define VER_2D

namespace SBV
{
#ifdef VER_2D
using Point = Eigen::Vector2d;
#else
using Point = Eigen::Vector3d;
#endif
    struct Mesh{
        std::vector<Point> vertices;
        std::vector<Eigen::Vector3i> triangles;
    };

    struct Curve{
        std::vector<Point> vertices;
        std::vector<Eigen::Vector2i> lines;
    };
}



#endif
