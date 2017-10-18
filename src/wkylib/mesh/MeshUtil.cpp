#include "MeshUtil.h"

using namespace zjucad::matrix;

namespace WKYLIB
{
    namespace Mesh
    {
        void computeNormal(const matrixr_t &vertices, const matrixs_t &triangles, matrixr_t &normals)
        {
            normals.resize(3, vertices.size(2));
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

        void computeNormal2D(const matrixr_t& vertices, const matrixr_t& lines, matrixr_t& normals)
        {
            if(vertices.size(1) != 2 || lines.size(1) != 2)
            {
                throw std::invalid_argument("You can only use 2d matrices for computeNormal2D");
            }
            normals.resize(2, vertices.size(2));
            matrixr_t vnSum = zeros(2, vertices.size(2));
            for(int i = 0; i < lines.size(2); i++){
                const matrixr_t& a = vertices(colon(), lines(0,i));
                const matrixr_t& b = vertices(colon(), lines(1,i));
                matrixr_t ab = b - a;

                matrixr_t faceNormal(2, 1);
                faceNormal[0] = ab[1];
                faceNormal[1] = -ab[0];
                for(int j = 0; j < 2; j++){
                    vnSum(colon(), lines(j,i)) += faceNormal;
                }
            }
            for(int i = 0; i < vertices.size(2); i++){
                double temp = norm(vnSum(colon(),i));
                normals(colon(), i) = vnSum(colon(),i) / temp;
            }
        }
    }
}
