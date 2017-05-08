#ifndef WKY_GEOMETRY_UTIL_H
#define WKY_GEOMETRY_UTIL_H

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int
#if (defined __GNUC__) && (__GNUC__>4 || __GNUC_MINOR__>=7)
  #undef _GLIBCXX_ATOMIC_BUILTINS
  #undef _GLIBCXX_USE_INT128
#endif

#include "Common.h"
#include <eigen3/Eigen/Dense>

namespace WKYLIB
{
    namespace Geometry
    {
        //zju matrix to eigen
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

        WKY_API void zju_mat_to_eigen(const matrixr_t& zju_matrix, Eigen::MatrixXd& output_eigen);

        //eigen to zju matrix
        WKY_API void eigen_to_zju_mat(const Eigen::MatrixXd& eigen, matrixr_t& output_matrix);

        WKY_API void eigen_to_zju_mat(const Eigen::MatrixXi& eigen, matrixs_t& output_matrix);

        WKY_API void eigen_to_zju_vector(const Eigen::VectorXd& eigen, matrixr_t& output_vector);

        WKY_API void eigen_to_zju_vector(const Eigen::VectorXi& eigen, matrixs_t& output_vector);

        WKY_API void eigen_to_zju_vector2d(const Eigen::Vector2d& eigen, matrixr_t& output_vector);

        WKY_API void eigen_to_zju_vector2d(const Eigen::Vector2i& eigen, matrixs_t& output_vector);

        WKY_API void eigen_to_zju_vector3d(const Eigen::Vector3d& eigen, matrixr_t& output_vector);

        WKY_API void eigen_to_zju_vector3d(const Eigen::Vector3i& eigen, matrixs_t& output_vector);

        WKY_API void eigen_to_zju_vector4d(const Eigen::Vector4d& eigen, matrixr_t& output_vector);

        WKY_API void eigen_to_zju_vector4d(const Eigen::Vector4i& eigen, matrixs_t& output_vector);
    }
}

#endif
