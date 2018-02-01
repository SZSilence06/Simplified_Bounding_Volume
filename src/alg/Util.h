#ifndef WKY_SBV_UTIL_H
#define WKY_SBV_UTIL_H

#include "Common.h"
#include <eigen3/Eigen/Dense>

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif

namespace SBV
{
    class Util{
    public:
        static double integrateOverTriangle(const vec3_t& x, const mat3x3_t &triangle);

        __host__ __device__ static double GPU_integrateOverTriangle(const Eigen::Vector3d& x, const Eigen::Matrix3d &triangle);
    };
}

#endif
