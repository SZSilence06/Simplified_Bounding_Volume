#include "MeshUtil.h"

using namespace zjucad::matrix;

namespace WKYLIB
{
    namespace Mesh
    {
        void computeNormal(const std::vector<Eigen::Vector3d>& vertices, const std::vector<Eigen::Vector3i>&triangles,
                           std::vector<Eigen::Vector3d>& normals)
        {
            std::vector<Eigen::Vector3d> vnSum(vertices.size());
            normals.resize(vertices.size());
            for(size_t i = 0; i < triangles.size(); i++){
                const Eigen::Vector3d& a = vertices[triangles[i][0]];
                const Eigen::Vector3d& b = vertices[triangles[i][1]];
                const Eigen::Vector3d& c = vertices[triangles[i][2]];

                Eigen::Vector3d faceNormal;
                faceNormal = (a-b).cross(a-c);
                for(int j = 0; j < 3; j++){
                    vnSum[triangles[i][j]] += faceNormal;
                }
            }
            for(size_t i = 0; i < vertices.size(); i++){
                double temp = vnSum[i].norm();
                normals[i] = vnSum[i] / temp;
            }
        }

        void computeNormal2D(const std::vector<Eigen::Vector2d>& vertices, const std::vector<Eigen::Vector2i>&lines,
                             std::vector<Eigen::Vector2d>& normals)
        {
            std::vector<Eigen::Vector2d> vnSum(vertices.size());
            normals.resize(vertices.size());
            for(size_t i = 0; i < lines.size(); i++){
                const Eigen::Vector2d& a = vertices[lines[i][0]];
                const Eigen::Vector2d& b = vertices[lines[i][1]];
                Eigen::Vector2d ab = b - a;

                Eigen::Vector2d faceNormal;
                faceNormal[0] = ab[1];
                faceNormal[1] = -ab[0];
                for(int j = 0; j < 2; j++){
                    vnSum[lines[i][j]] += faceNormal;
                }
            }
            for(size_t i = 0; i < vertices.size(); i++){
                double temp = vnSum[i].norm();
                normals[i] = vnSum[i] / temp;
            }
        }
    }
}
