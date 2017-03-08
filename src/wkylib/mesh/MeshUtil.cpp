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
                matrixr_t a = vertices(colon(),triangles(0,i));
                matrixr_t b = vertices(colon(),triangles(1,i));
                matrixr_t c = vertices(colon(),triangles(2,i));

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

        void computeCurveNormal2D(const matrixr_t& vertices, matrixr_t& normals)
        {
            if(vertices.size(2) < 3 || vertices.size(1) != 2)
            {
                throw std::invalid_argument("computeCurveNormal2D arg invalid.");
            }

            matrixr_t vnSum = zeros(2, vertices.size(2));
            for(int i = 0; i < vertices.size(2); i++){
                matrixr_t a;
                matrixr_t b;

                if(i < vertices.size(2) - 1)
                {
                    a = vertices(colon(), i);
                    b = vertices(colon(), i + 1);
                }
                else
                {
                    a = vertices(colon(), i);
                    b = vertices(colon(), 0);
                }

                matrixr_t ab = b - a;
                matrixr_t faceNormal(2, 1);
                faceNormal[0] = -ab[1];
                faceNormal[1] = ab[0];

                if(i < vertices.size(2) - 1)
                {
                    vnSum(colon(), i) += faceNormal;
                    vnSum(colon(), i + 1) += faceNormal;
                }
                else
                {
                    vnSum(colon(), i) += faceNormal;
                    vnSum(colon(), 0) += faceNormal;
                }
            }
            for(int i = 0; i < vertices.size(2); i++){
                double temp = norm(vnSum(colon(),i));
                normals(colon(), i) = vnSum(colon(),i) / temp;
            }
        }
    }
}
