#include "Simplifier.h"
#include "Refinement.h"
#include <wkylib/mesh/MeshUtil.h>
#include <wkylib/mesh/IO.h>
#include <jtflib/mesh/io.h>

using namespace zjucad::matrix;

namespace SBV
{
    Simplifier::Simplifier(Curve &mesh) : mSourceMesh(mesh)
    {

    }

    void Simplifier::simplify()
    {
        generateShells();
        //refine();
    }

    void Simplifier::generateShells()
    {
        matrixr_t normals(2, mSourceMesh.vertices.size(2));

        WKYLIB::Mesh::computeNormal2D(mSourceMesh.vertices, mSourceMesh.lines, normals);

        matrixr_t shell(2, mSourceMesh.vertices.size(2));

        for(int i = 0; i < mSourceMesh.vertices.size(2); i++)
        {
            const matrixr_t& vertex = mSourceMesh.vertices(colon(), i);
            shell(colon(), i) = vertex + normals(colon(), i) * mMaxDistance;
        }

        std::vector<matrixr_t> sampled_shell;

        sample(shell, mSourceMesh.lines, sampled_shell);
        mOuterShell.resize(2, sampled_shell.size());
        for(int i = 0; i < mOuterShell.size(2); i++)
        {
            mOuterShell(colon(), i) = sampled_shell[i];
        }

        sample(mSourceMesh.vertices, mSourceMesh.lines, sampled_shell);
        mInnerShell.resize(2, sampled_shell.size());
        for(int i = 0; i < mInnerShell.size(2); i++)
        {
            mInnerShell(colon(), i) = sampled_shell[i];
        }

        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writePoints2D(mInnerShell, mOutputDirectory + "/inner_shell.vtk");
            WKYLIB::Mesh::writePoints2D(mOuterShell, mOutputDirectory + "/outer_shell.vtk");
        }
    }

    void Simplifier::sample(const matrixr_t &vertices, const matrixs_t &triangles, std::vector<matrixr_t> &output_samples)
    {
        output_samples.clear();
        for(int i = 0; i < triangles.size(2); i++)
        {
            const matrixr_t& a = vertices(colon(), triangles(0, i));
            const matrixr_t& b = vertices(colon(), triangles(1, i));

            matrixr_t ab = b - a;

            int sample_count_ab = (int)(norm(ab) / mSampleRadius);

            ab /= norm(ab);

            for(int i = 0; i < sample_count_ab; i++)
            {
                output_samples.push_back(a + mSampleRadius * i * ab);
            }
        }
    }

    void Simplifier::refine()
    {
        Refinement refinement(mInnerShell, mOuterShell);
        std::vector<size_t> refineOutput;
        refinement.refine(refineOutput);
    }
}
