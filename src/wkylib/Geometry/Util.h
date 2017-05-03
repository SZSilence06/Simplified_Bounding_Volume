#ifndef WKY_GEOMETRY_UTIL_H
#define WKY_GEOMETRY_UTIL_H

#include "Common.h"
#include <eigen3/Eigen/Dense>

namespace WKYLIB
{
    namespace Geometry
    {
        WKY_API void zju_mat_to_eigen(const matrixr_t& zju_matrix, Eigen::MatrixXd& output_eigen);

        WKY_API void zju_mat_to_eigen(const matrixs_t& zju_matrix, Eigen::MatrixXi& output_eigen);

        WKY_API void zju_vector_to_eigen(const matrixr_t& zju_vector, Eigen::VectorXd& output_eigen);

        WKY_API void zju_vector_to_eigen(const matrixs_t& zju_vector, Eigen::VectorXi& output_eigen);

        WKY_API void zju_vector2d_to_eigen(const matrixr_t& zju_vector, Eigen::Vector2d& output_eigen);

        WKY_API void zju_vector2d_to_eigen(const matrixs_t& zju_vector, Eigen::Vector2i& output_eigen);

        WKY_API void zju_vector3d_to_eigen(const matrixr_t& zju_vector, Eigen::Vector3d& output_eigen);

        WKY_API void zju_vector3d_to_eigen(const matrixs_t& zju_vector, Eigen::Vector3i& output_eigen);

        WKY_API void zju_vector4d_to_eigen(const matrixr_t& zju_vector, Eigen::Vector4d& output_eigen);

        WKY_API void zju_vector4d_to_eigen(const matrixs_t& zju_vector, Eigen::Vector4i& output_eigen);
    }
}

#endif
