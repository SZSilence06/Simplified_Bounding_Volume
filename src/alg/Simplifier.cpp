#include "Simplifier.h"
#include "Refinement.h"
#include "BoundaryCollapse.h"
#include <wkylib/mesh/MeshUtil.h>
#include <wkylib/mesh/IO.h>
#include <wkylib/debug_util.h>
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

        std::cout << "Start refinement..." << std::endl;
        refine();

        std::cout << "Start Boundary collapse..." << std::endl;
        collapseBoundary();

        std::cout << "Start mutual tessellation..." << std::endl;
        mutualTessellate();
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
            WKYLIB::Mesh::writePoints2D(mOutputDirectory + "/inner_shell.vtk", mInnerShell);
            WKYLIB::Mesh::writePoints2D(mOutputDirectory + "/outer_shell.vtk", mOuterShell);
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
        WKYLIB::DebugTimer timer("refinement");
        timer.start();
        Refinement refinement(mInnerShell, mOuterShell, mTriangulation, mAlpha, mSampleRadius);
        refinement.refine();
        timer.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/refined_shell.obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/refined_zero_set.obj", mTriangulation.getZeroSet().vertices,
                                       mTriangulation.getZeroSet().lines);
        }
    }

    void Simplifier::collapseBoundary()
    {
        WKYLIB::DebugTimer timer("Boundary Collapse");
        timer.start();
        BoundaryCollapse collapser(mTriangulation, mInnerShell, mOuterShell);
        collapser.collapse();
        timer.end();

         mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/boundary_collapsed_shell.obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/boundary_collapsed_zero_set.obj", mTriangulation.getZeroSet().vertices,
                                       mTriangulation.getZeroSet().lines);
        }
    }

    void Simplifier::mutualTessellate()
    {
        WKYLIB::DebugTimer timer("Mutual Tesselation");
        timer.start();
        mTriangulation.mutualTessellate();
        timer.end();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/mutual_tesellation.obj", mTriangulation.vertices, mTriangulation.triangles);
        }
    }
}
