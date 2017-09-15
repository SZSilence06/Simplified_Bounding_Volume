#include "Sampler.h"
#include <wkylib/mesh/MeshUtil.h>
#include "External/poisson_disk/poisson_disk_wrapper/utils_sampling.hpp"

using namespace Utils_sampling;

namespace SBV
{
    void Sampler::poissonDisk(const matrixr_t &vertices, const matrixs_t &triangles, double radius, matrixr_t &samples, matrixr_t& sample_normals)
    {
        std::vector<Vec3> verts;
        for(int i = 0; i < vertices.size(2); i++)
        {
            Vec3 vert;
            vert.x = vertices(0, i);
            vert.y = vertices(1, i);
            vert.z = vertices(2, i);
            verts.push_back(vert);
        }
        std::vector<int> tris;
        for(int i = 0; i < triangles.size(2); i++)
        {
            for(int j = 0; j < 3; j++)
            {
                tris.push_back(triangles(j, i));
            }
        }

        matrixr_t normals;
        WKYLIB::Mesh::computeNormal(vertices, triangles, normals);
        std::vector<Vec3> nors;
        for(int i = 0; i < normals.size(2); i++)
        {
            Vec3 nor;
            nor.x = normals(0, i);
            nor.y = normals(1, i);
            nor.z = normals(2, i);
            nors.push_back(nor);
        }

        std::vector<Vec3> out_samples;
        std::vector<Vec3> out_sample_normals;
        poisson_disk(radius, 0, verts, nors, tris, out_samples, out_sample_normals);

        samples.resize(3, out_samples.size());
        for(int i = 0; i < out_samples.size(); i++)
        {
            samples(0, i) = out_samples[i].x;
            samples(1, i) = out_samples[i].y;
            samples(2, i) = out_samples[i].z;
        }

        sample_normals.resize(3, out_sample_normals.size());
        for(int i = 0; i < out_samples.size(); i++)
        {
            sample_normals(0, i) = out_sample_normals[i].x;
            sample_normals(1, i) = out_sample_normals[i].y;
            sample_normals(2, i) = out_sample_normals[i].z;
        }
    }
}
