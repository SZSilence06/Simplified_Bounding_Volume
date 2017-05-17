#include "MeshUtil.h"

using namespace zjucad::matrix;

namespace WKYLIB
{
    namespace Mesh
    {
        void computeNormal(const matrixr_t &vertices, const matrixs_t &triangles, matrixr_t &normals)
        {
            matrixr_t vnSum = zeros(3, vertices.size(2));
            for(int i = 0; i < triangles.size(2); i++){
                const matrixr_t& a = vertices(colon(),triangles(0,i));
                const matrixr_t& b = vertices(colon(),triangles(1,i));
                const matrixr_t& c = vertices(colon(),triangles(2,i));

                matrixr_t faceNormal(3, 1);
                faceNormal = cross(a-b,a-c);
                for(int j = 0; j < 3; j++){
                    vnSum(colon(), triangles(j,i)) += faceNormal;
                }
            }
            for(int i = 0; i < vertices.size(2); i++){
                double temp = norm(vnSum(colon(),i));
                normals(colon(), i) = vnSum(colon(),i) / temp;
            }
        }

        void computeNormal2D(const std::vector<Eigen::Vector2d>& vertices, const std::vector<Eigen::Vector2i>&lines,
                             std::vector<Eigen::Vector2d>& normals)
        {
            std::vector<Eigen::Vector2d> vnSum(vertices.size());
            normals.resize(vertices.size());
            for(int i = 0; i < lines.size(); i++){
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
            for(int i = 0; i < vertices.size(); i++){
                double temp = vnSum[i].norm();
                normals[i] = vnSum[i] / temp;
            }
        }
    }
}
