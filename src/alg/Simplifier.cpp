#include "Simplifier.h"
#include <wkylib/mesh/MeshUtil.h>
#include <wkylib/mesh/IO.h>
#include <jtflib/mesh/io.h>

using namespace zjucad::matrix;

namespace SBV
{
    Simplifier::Simplifier(Mesh &mesh) : mSourceMesh(mesh)
    {

    }

    void Simplifier::simplify()
    {
        generateShells();
    }

    void Simplifier::generateShells()
    {
        matrixr_t normals(3, mSourceMesh.vertices.size(2));

        WKYLIB::Mesh::computeNormal(mSourceMesh.vertices, mSourceMesh.triangles, normals);

        matrixr_t shell(3, mSourceMesh.vertices.size(2));

        for(int i = 0; i < mSourceMesh.vertices.size(2); i++)
        {
            matrixr_t vertex = mSourceMesh.vertices(colon(), i);
            shell(colon(), i) = vertex + normals(colon(), i) * mMaxDistance;
        }

        std::vector<matrixr_t> sampled_shell;

        sample(shell, mSourceMesh.triangles, sampled_shell);
        mOuterShell.resize(3, sampled_shell.size());
        for(int i = 0; i < mOuterShell.size(2); i++)
        {
            mOuterShell(colon(), i) = sampled_shell[i];
        }

        sample(mSourceMesh.vertices, mSourceMesh.triangles, sampled_shell);
        mInnerShell.resize(3, sampled_shell.size());
        for(int i = 0; i < mInnerShell.size(2); i++)
        {
            mInnerShell(colon(), i) = sampled_shell[i];
        }

        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writePoints(mInnerShell, mOutputDirectory + "/inner_shell.vtk");
            WKYLIB::Mesh::writePoints(mOuterShell, mOutputDirectory + "/outer_shell.vtk");
        }
    }

    void Simplifier::sample(const matrixr_t &vertices, const matrixs_t &triangles, std::vector<matrixr_t> &output_samples)
    {
        output_samples.clear();
        for(int i = 0; i < triangles.size(2); i++)
        {
            const matrixr_t& a = vertices(colon(), triangles(0, i));
            const matrixr_t& b = vertices(colon(), triangles(1, i));
            const matrixr_t& c = vertices(colon(), triangles(2, i));

            matrixr_t ab = b - a;
            matrixr_t ac = c - a;

            int sample_count_ab = (int)(norm(ab) / mSampleRadius);
            int sample_count_ac = (int)(norm(ac) / mSampleRadius);

            ab /= norm(ab);
            ac /= norm(ac);

            for(int i = 0; i < sample_count_ab; i++)
            {
                for(int j = 0; j < sample_count_ac; j++)
                {
                    output_samples.push_back(a + mSampleRadius * i * ab + mSampleRadius * j * ac);
                }
            }
        }
    }
}
