#include "Simplifier.h"
#include "Refinement.h"
#include "EdgeCollapse.h"
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

        WKYLIB::DebugTimer timer("Build CudaController");
        timer.start();
        mCudaController.build(mShell.mInnerShell, mShell.mOuterShell, mTriangulation);
        timer.end();

        std::cout << "Start Boundary collapse..." << std::endl;
        collapseBoundary();

        std::cout << "Start mutual tessellation..." << std::endl;
        mutualTessellate();

        std::cout << "Start zero set collapse..." << std::endl;
        collapseZeroSet();
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
        mShell.mOuterShell.resize(2, sampled_shell.size());
        for(int i = 0; i < mShell.mOuterShell.size(2); i++)
        {
            mShell.mOuterShell(colon(), i) = sampled_shell[i];
        }

        sample(mSourceMesh.vertices, mSourceMesh.lines, sampled_shell);
        mShell.mInnerShell.resize(2, sampled_shell.size());
        for(int i = 0; i < mShell.mInnerShell.size(2); i++)
        {
            mShell.mInnerShell(colon(), i) = sampled_shell[i];
        }

        mShell.buildKdTree();

        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writePoints2D(mOutputDirectory + "/inner_shell.vtk", mShell.mInnerShell);
            WKYLIB::Mesh::writePoints2D(mOutputDirectory + "/outer_shell.vtk", mShell.mOuterShell);
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
        Refinement refinement(mShell, mTriangulation, mAlpha, mSampleRadius);
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
        WKYLIB::DebugTimer timer("Boundary Collapse(Half Edge)");
        timer.start();
        EdgeCollapse collapserHalfEdge(mTriangulation, mShell, mCudaController, EdgeCollapse::BOUNDARY, true, mSampleRadius);
        collapserHalfEdge.collapse();
        timer.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/boundary_collapsed_shell(half_edge).obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/boundary_collapsed_zero_set(half_edge).obj", mTriangulation.getZeroSet().vertices,
                                       mTriangulation.getZeroSet().lines);
        }

        WKYLIB::DebugTimer timer2("Boundary Collapse(General)");
        timer2.start();
        EdgeCollapse collapserGeneral(mTriangulation, mShell, mCudaController, EdgeCollapse::BOUNDARY, false, mSampleRadius);
        collapserGeneral.collapse();
        timer2.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/boundary_collapsed_shell(general).obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/boundary_collapsed_zero_set(general).obj", mTriangulation.getZeroSet().vertices,
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

    void Simplifier::collapseZeroSet()
    {
        WKYLIB::DebugTimer timer("Zero Set Collapse(Half Edge)");
        timer.start();
        EdgeCollapse collapserHalfEdge(mTriangulation, mShell, mCudaController, EdgeCollapse::ZERO_SET, true, mSampleRadius);
        collapserHalfEdge.collapse();
        timer.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/zero_set_collapsed_shell.obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/zero_set_collapsed_zero_set.obj", mTriangulation.getZeroSet().vertices,
                                       mTriangulation.getZeroSet().lines);
        }

        WKYLIB::DebugTimer timer2("Zero Set Collapse(General)");
        timer2.start();
        EdgeCollapse collapserGeneral(mTriangulation, mShell, mCudaController, EdgeCollapse::ZERO_SET, false, mSampleRadius);
        collapserGeneral.collapse();
        timer2.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/zero_set_collapsed_shell(general).obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/zero_set_collapsed_zero_set(general).obj", mTriangulation.getZeroSet().vertices,
                                       mTriangulation.getZeroSet().lines);
        }
    }
}
