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
        mTimerSimplify.start();

        generateShells();

        std::cout << "Start refinement..." << std::endl;
        refine();

        std::cout << "Start Boundary collapse..." << std::endl;
        collapseBoundary();

        std::cout << "Start mutual tessellation..." << std::endl;
        mutualTessellate();

        std::cout << "Start zero set collapse..." << std::endl;
        collapseZeroSet();

        mTimerSimplify.end();

        if(mNeedGenTempResult)
            writeSummary();
    }

    void Simplifier::generateShells()
    {
        std::vector<Point> normals;

        WKYLIB::Mesh::computeNormal2D(mSourceMesh.vertices, mSourceMesh.lines, normals);

        std::vector<Point> shell;
        shell.reserve(mSourceMesh.vertices.size());

        for(size_t i = 0; i < mSourceMesh.vertices.size(); i++)
        {
            const Point& vertex = mSourceMesh.vertices[i];
            shell.push_back(vertex + normals[i] * mMaxDistance);
        }

        std::vector<Point> sampled_shell;

        sample(shell, mSourceMesh.lines, sampled_shell);
        mShell.mOuterShell.reserve(sampled_shell.size());
        for(size_t i = 0; i < sampled_shell.size(); i++)
        {
            mShell.mOuterShell.push_back(sampled_shell[i]);
        }

        sample(mSourceMesh.vertices, mSourceMesh.lines, sampled_shell);
        mShell.mInnerShell.reserve(sampled_shell.size());
        for(size_t i = 0; i < sampled_shell.size(); i++)
        {
            mShell.mInnerShell.push_back(sampled_shell[i]);
        }

        mShell.buildKdTree();

        if(mNeedGenTempResult)
        {
#ifdef VER_2D
            WKYLIB::Mesh::writePoints2D(mOutputDirectory + "/inner_shell.vtk", mShell.mInnerShell);
            WKYLIB::Mesh::writePoints2D(mOutputDirectory + "/outer_shell.vtk", mShell.mOuterShell);
#else
            WKYLIB::Mesh::writePoints(mOutputDirectory + "/inner_shell.vtk", mShell.mInnerShell);
            WKYLIB::Mesh::writePoints(mOutputDirectory + "/outer_shell.vtk", mShell.mOuterShell);
#endif
        }
    }

#ifdef VER_2D
    void Simplifier::sample(const std::vector<Point> &vertices, const std::vector<Eigen::Vector2i> &lines,
                            std::vector<Point> &output_samples)
    {
        output_samples.clear();
        for(size_t i = 0; i < lines.size(); i++)
        {
            const Point& a = vertices[lines[i][0]];
            const Point& b = vertices[lines[i][1]];

            Point ab = b - a;

            int sample_count_ab = (int)(ab.norm() / mSampleRadius);

            ab /= ab.norm();

            for(int i = 0; i < sample_count_ab; i++)
            {
                output_samples.push_back(a + mSampleRadius * i * ab);
            }
        }
    }
#endif

    void Simplifier::refine()
    {
        mTimerRefine.start();
        Refinement refinement(mShell, mTriangulation, mAlpha, mSampleRadius);
        refinement.refine();
        mTimerRefine.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
#ifdef VER_2D
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/refined_shell.obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/refined_zero_set.obj", mTriangulation.getZeroSet().vertices,
                                       mTriangulation.getZeroSet().lines);
#else
#endif
        }
    }

    void Simplifier::collapseBoundary()
    {
        mTimerBoundaryHalfEdge.start();
        EdgeCollapse collapserHalfEdge(mTriangulation, mShell, mCudaController, EdgeCollapse::BOUNDARY, true, mSampleRadius);
        collapserHalfEdge.collapse();
        mTimerBoundaryHalfEdge.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
#ifdef VER_2D
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/boundary_collapsed_shell(half_edge).obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/boundary_collapsed_zero_set(half_edge).obj", mTriangulation.getZeroSet().vertices,
                                       mTriangulation.getZeroSet().lines);
#else
#endif
        }

        mTimerBuildCuda.start();
        mCudaController.build(mShell.mInnerShell, mShell.mOuterShell, mTriangulation);
        mTimerBuildCuda.suspend();

        mTimerBoundaryGeneral.start();
        EdgeCollapse collapserGeneral(mTriangulation, mShell, mCudaController, EdgeCollapse::BOUNDARY, false, mSampleRadius);
        collapserGeneral.collapse();
        mTimerBoundaryGeneral.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
#ifdef VER_2D
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/boundary_collapsed_shell(general).obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/boundary_collapsed_zero_set(general).obj", mTriangulation.getZeroSet().vertices,
                                       mTriangulation.getZeroSet().lines);
#else
#endif
        }
    }

    void Simplifier::mutualTessellate()
    {
        mTimerMutualTessellation.start();
        mTriangulation.mutualTessellate();
        mTimerMutualTessellation.end();

        if(mNeedGenTempResult)
        {
#ifdef VER_2D
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/mutual_tesellation.obj", mTriangulation.vertices, mTriangulation.triangles);
#else
#endif
        }
    }

    void Simplifier::collapseZeroSet()
    {
        mTimerZeroSetHalfEdge.start();
        EdgeCollapse collapserHalfEdge(mTriangulation, mShell, mCudaController, EdgeCollapse::ZERO_SET, true, mSampleRadius);
        collapserHalfEdge.collapse();
        mTimerZeroSetHalfEdge.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
#ifdef VER_2D
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/zero_set_collapsed_shell.obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/zero_set_collapsed_zero_set.obj", mTriangulation.getZeroSet().vertices,
                                       mTriangulation.getZeroSet().lines);
#else
#endif
        }

        mTimerBuildCuda.resume();
        mCudaController.build(mShell.mInnerShell, mShell.mOuterShell, mTriangulation);
        mTimerBuildCuda.end();

        mTimerZeroSetGeneral.start();
        EdgeCollapse collapserGeneral(mTriangulation, mShell, mCudaController, EdgeCollapse::ZERO_SET, false, mSampleRadius);
        collapserGeneral.collapse();
        mTimerZeroSetGeneral.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
#ifdef VER_2D
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/zero_set_collapsed_shell(general).obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/zero_set_collapsed_zero_set(general).obj", mTriangulation.getZeroSet().vertices,
                                       mTriangulation.getZeroSet().lines);
#else
#endif
        }
    }

    void Simplifier::writeSummary()
    {
        logTimer("Refinement : ", mTimerRefine);
        logTimer("Build CudaController : ", mTimerBuildCuda);
        logTimer("Boundary Collapse(Half Edge) : ", mTimerBoundaryHalfEdge);
        logTimer("Boundary Collapse(General) : ", mTimerBoundaryGeneral);
        logTimer("Mutual Tessellation : ", mTimerMutualTessellation);
        logTimer("Zero Set Collapse(Half Edge) : ", mTimerZeroSetHalfEdge);
        logTimer("Zero Set Collapse(General) : ", mTimerZeroSetGeneral);
        logTimer("Total : ", mTimerSimplify);
    }

    void Simplifier::logTimer(const std::string& prefix, const DebugTimer &timer)
    {
        Logger& logger = Logger::getInstance();
        logger.log(prefix + std::to_string(timer.getTime()) + " ms.");
    }
}
