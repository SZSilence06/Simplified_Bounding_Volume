#include "Util.h"

namespace WKYLIB
{
    namespace Geometry
    {
        //zju matrix to eigen
        WKY_API void zju_mat_to_eigen(const matrixr_t &zju_matrix, Eigen::MatrixXd &output_eigen)
        {
            output_eigen.resize(zju_matrix.size(1), zju_matrix.size(2));
            for(int i = 0; i < zju_matrix.size(1); i++)
            {
                for(int j = 0; j < zju_matrix.size(2); j++)
                {
                    output_eigen(i, j) = zju_matrix(i, j);
                }
            }
        }

        WKY_API void zju_mat_to_eigen(const matrixs_t &zju_matrix, Eigen::MatrixXi &output_eigen)
        {
            output_eigen.resize(zju_matrix.size(1), zju_matrix.size(2));
            for(int i = 0; i < zju_matrix.size(1); i++)
            {
                for(int j = 0; j < zju_matrix.size(2); j++)
                {
                    output_eigen(i, j) = zju_matrix(i, j);
                }
            }
        }

        WKY_API void zju_vector_to_eigen(const matrixr_t &zju_vector, Eigen::VectorXd &output_eigen)
        {
            if(zju_vector.size(2) != 1)
            {
                throw std::invalid_argument("you can only use vectors in zju_vector_to_eigen().");
            }

            output_eigen.resize(1, zju_vector.size(2));
            for(int i = 0; i < zju_vector.size(2); i++)
            {
                output_eigen[i] = zju_vector[i];
            }
        }

        WKY_API void zju_vector_to_eigen(const matrixs_t &zju_vector, Eigen::VectorXi &output_eigen)
        {
            if(zju_vector.size(2) != 1)
            {
                throw std::invalid_argument("you can only use vectors in zju_vector_to_eigen().");
            }

            output_eigen.resize(zju_vector.size(), 1);
            for(int i = 0; i < zju_vector.size(); i++)
            {
                output_eigen[i] = zju_vector[i];
            }
        }

        WKY_API void zju_vector2d_to_eigen(const matrixr_t &zju_vector, Eigen::Vector2d &output_eigen)
        {
            if(zju_vector.size(2) != 1 || zju_vector.size(1) != 2)
            {
                throw std::invalid_argument("you can only use 2d vectors in zju_vector2d_to_eigen().");
            }

            output_eigen.resize(2, 1);
            for(int i = 0; i < 2; i++)
            {
                output_eigen[i] = zju_vector[i];
            }
        }

        WKY_API void zju_vector2d_to_eigen(const matrixs_t &zju_vector, Eigen::Vector2i &output_eigen)
        {
            if(zju_vector.size(2) != 1 || zju_vector.size(1) != 2)
            {
                throw std::invalid_argument("you can only use 2d vectors in zju_vector2d_to_eigen().");
            }

            output_eigen.resize(2, 1);
            for(int i = 0; i < 2; i++)
            {
                output_eigen[i] = zju_vector[i];
            }
        }

        WKY_API void zju_vector3d_to_eigen(const matrixr_t &zju_vector, Eigen::Vector3d &output_eigen)
        {
            if(zju_vector.size(2) != 1 || zju_vector.size(1) != 3)
            {
                throw std::invalid_argument("you can only use 3d vectors in zju_vector3d_to_eigen().");
            }
            output_eigen.resize(3, 1);
            for(int i = 0; i < 3; i++)
            {
                output_eigen[i] = zju_vector[i];
            }
        }

        WKY_API void zju_vector3d_to_eigen(const matrixs_t &zju_vector, Eigen::Vector3i &output_eigen)
        {
            if(zju_vector.size(2) != 1 || zju_vector.size(1) != 3)
            {
                throw std::invalid_argument("you can only use 3d vectors in zju_vector3d_to_eigen().");
            }

            output_eigen.resize(3, 1);
            for(int i = 0; i < 3; i++)
            {
                output_eigen[i] = zju_vector[i];
            }
        }

        WKY_API void zju_vector4d_to_eigen(const matrixr_t &zju_vector, Eigen::Vector4d &output_eigen)
        {
            if(zju_vector.size(2) != 1 || zju_vector.size(1) != 4)
            {
                throw std::invalid_argument("you can only use 4d vectors in zju_vector4d_to_eigen().");
            }

            output_eigen.resize(4, 1);
            for(int i = 0; i < 4; i++)
            {
                output_eigen[i] = zju_vector[i];
            }
        }

        WKY_API void zju_vector4d_to_eigen(const matrixs_t &zju_vector, Eigen::Vector4i &output_eigen)
        {
            if(zju_vector.size(2) != 1 || zju_vector.size(1) != 4)
            {
                throw std::invalid_argument("you can only use 4d vectors in zju_vector4d_to_eigen().");
            }

            output_eigen.resize(4, 1);
            for(int i = 0; i < 4; i++)
            {
                output_eigen[i] = zju_vector[i];
            }
        }

        //eigen to zju matrix
        WKY_API void eigen_to_zju_mat(const Eigen::MatrixXd& eigen, matrixr_t& output_matrix)
        {
            output_matrix.resize(eigen.rows(), eigen.cols());
            for(int i = 0; i < eigen.rows(); i++)
            {
                for(int j = 0; j < eigen.cols(); j++)
                {
                    output_matrix(i, j) = eigen(i, j);
                }
            }
        }

        WKY_API void eigen_to_zju_mat(const Eigen::MatrixXi& eigen, matrixs_t& output_matrix)
        {
            output_matrix.resize(eigen.rows(), eigen.cols());
            for(int i = 0; i < eigen.rows(); i++)
            {
                for(int j = 0; j < eigen.cols(); j++)
                {
                    output_matrix(i, j) = eigen(i, j);
                }
            }
        }

        WKY_API void eigen_to_zju_vector(const Eigen::VectorXd& eigen, matrixr_t& output_vector)
        {
            output_vector.resize(eigen.rows(), 1);
            for(int i = 0; i < eigen.rows(); i++)
            {
                output_vector[i] = eigen[i];
            }
        }

        WKY_API void eigen_to_zju_vector(const Eigen::VectorXi& eigen, matrixs_t& output_vector)
        {
            output_vector.resize(eigen.rows(), 1);
            for(int i = 0; i < eigen.rows(); i++)
            {
                output_vector[i] = eigen[i];
            }
        }

        WKY_API void eigen_to_zju_vector2d(const Eigen::Vector2d& eigen, matrixr_t& output_vector)
        {
            output_vector.resize(eigen.rows(), 1);
            for(int i = 0; i < eigen.rows(); i++)
            {
                output_vector[i] = eigen[i];
            }
        }

        WKY_API void eigen_to_zju_vector2d(const Eigen::Vector2i& eigen, matrixs_t& output_vector)
        {
            output_vector.resize(eigen.rows(), 1);
            for(int i = 0; i < eigen.rows(); i++)
            {
                output_vector[i] = eigen[i];
            }
        }

        WKY_API void eigen_to_zju_vector3d(const Eigen::Vector3d& eigen, matrixr_t& output_vector)
        {
            output_vector.resize(eigen.rows(), 1);
            for(int i = 0; i < eigen.rows(); i++)
            {
                output_vector[i] = eigen[i];
            }
        }

        WKY_API void eigen_to_zju_vector3d(const Eigen::Vector3i& eigen, matrixs_t& output_vector)
        {
            output_vector.resize(eigen.rows(), 1);
            for(int i = 0; i < eigen.rows(); i++)
            {
                output_vector[i] = eigen[i];
            }
        }

        WKY_API void eigen_to_zju_vector42d(const Eigen::Vector4d& eigen, matrixr_t& output_vector)
        {
            output_vector.resize(eigen.rows(), 1);
            for(int i = 0; i < eigen.rows(); i++)
            {
                output_vector[i] = eigen[i];
            }
        }

        WKY_API void eigen_to_zju_vector4d(const Eigen::Vector4i& eigen, matrixs_t& output_vector)
        {
            output_vector.resize(eigen.rows(), 1);
            for(int i = 0; i < eigen.rows(); i++)
            {
                output_vector[i] = eigen[i];
            }
        }
    }
}
