#ifndef WKY_SBV_COMMON_H
#define WKY_SBV_COMMON_H

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int

#include <vector>
#include "eigen3.3/Eigen/Dense"

#define VER_2D

#pragma NVCC diagnostic ignored "-Wall"
#if (defined __GNUC__) && (__GNUC__>4 || __GNUC_MINOR__>=7)
  #undef _GLIBCXX_ATOMIC_BUILTINS
  #undef _GLIBCXX_USE_INT128
#endif

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
