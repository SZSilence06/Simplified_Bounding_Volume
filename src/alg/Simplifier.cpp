#include "Simplifier.h"
#include "ShellGenerator.h"
#include "Refinement.h"
#include "EdgeCollapse.h"
#include "Sampler.h"
#include <wkylib/mesh/MeshUtil.h>
#include <wkylib/mesh/IO.h>
#include <wkylib/debug_util.h>
#include <jtflib/mesh/io.h>

using namespace zjucad::matrix;

namespace SBV
{
    void uniformSpherical(vec3_t& p)
    {
        while(true)
        {
            double x1 = rand() * (2.0 / RAND_MAX) - 1;
            double x2 = rand() * (2.0 / RAND_MAX) - 1;
            double temp = x1*x1 + x2*x2;
            if(temp < 1)
            {
                p[0] = 2 * x1 * sqrt(1 - temp);
                p[1] = 2 * x2 * sqrt(1 - temp);
                p[2] = 1 - 2 * temp;
                break;
            }
        }
    }

    Simplifier::Simplifier(Shell &shell) : mShell(shell)
    {
        srand(time(0));
    }

    void Simplifier::simplify()
    {
        mTimerSimplify.start();

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

    void Simplifier::sample(const matrixr_t &vertices, const matrixs_t &triangles, matrixr_t &output_samples, matrixr_t& output_normals)
    {
        Sampler::poissonDisk(vertices, triangles, mSampleRadius, output_samples, output_normals);
        /*output_samples.clear();
        for(int i = 0; i < triangles.size(2); i++)
        {
            const matrixr_t& a = vertices(colon(), triangles(0, i));
            const matrixr_t& b = vertices(colon(), triangles(1, i));
            const matrixr_t& c = vertices(colon(), triangles(2, i));

            matrixr_t ab = b - a;
            matrixr_t ac = c - a;
            double normAB = norm(ab);
            double normAC = norm(ac);

            ab /= normAB;
            ac /= normAC;

            double sample_count_ab = normAB / mSampleRadius;
            double sample_count_ac = normAC / mSampleRadius;

            for(int i = 1; i <= sample_count_ab; i++)
            {
                double baryAB = i / sample_count_ab;
                int max_sample_ac = (1 - baryAB) * sample_count_ac;
                for(int j = 1; j < max_sample_ac; j++)
                {
                    output_samples.push_back(a + mSampleRadius * i * ab + mSampleRadius * j * ac);
                }
            }
        }*/
    }

    void Simplifier::refine()
    {
        mTimerRefine.start();
        Refinement refinement(mShell, mTriangulation, mAlpha, mSampleRadius);
        refinement.refine();
        mTimerRefine.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeTetra(mOutputDirectory + "/refined_shell.vtk", mTriangulation.vertices, mTriangulation.cells);
            jtf::mesh::save_obj((mOutputDirectory + "/refined_zero_set.obj").c_str(), mTriangulation.getZeroSet().triangles,
                                       mTriangulation.getZeroSet().vertices);
            mTriangulation.outputBoundaryFaces(mOutputDirectory + "/refined");
        }
    }

    void Simplifier::collapseBoundary()
    {
        mTimerBoundaryHalfEdge.start();
        EdgeCollapse collapserHalfEdge(mTriangulation, mShell, EdgeCollapse::BOUNDARY, true, 1e-6);
        collapserHalfEdge.collapse();
        mTimerBoundaryHalfEdge.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeTetra(mOutputDirectory + "/boundary_collapsed_shell(half_edge).vtk", mTriangulation.vertices, mTriangulation.cells);
            jtf::mesh::save_obj((mOutputDirectory + "/boundary_collapsed_zero_set(half_edge).obj").c_str(), mTriangulation.getZeroSet().triangles,
                                       mTriangulation.getZeroSet().vertices);
            mTriangulation.outputBoundaryFaces(mOutputDirectory + "/boundary_collapsed(half_edge)");
        }

        mTimerBoundaryGeneral.start();
        EdgeCollapse collapserGeneral(mTriangulation, mShell, EdgeCollapse::BOUNDARY, false, 1e-6);
        collapserGeneral.collapse();
        mTimerBoundaryGeneral.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeTetra(mOutputDirectory + "/boundary_collapsed_shell(general).vtk", mTriangulation.vertices, mTriangulation.cells);
            jtf::mesh::save_obj((mOutputDirectory + "/boundary_collapsed_zero_set(general).obj").c_str(), mTriangulation.getZeroSet().triangles,
                                       mTriangulation.getZeroSet().vertices);
            mTriangulation.outputBoundaryFaces(mOutputDirectory + "/boundary_collapsed(general)");
        }
    }

    void Simplifier::mutualTessellate()
    {
        mTimerMutualTessellation.start();
        mTriangulation.mutualTessellate();
        mTimerMutualTessellation.end();

        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeTetra(mOutputDirectory + "/mutual_tessellation.vtk", mTriangulation.vertices, mTriangulation.cells);
        }
    }

    void Simplifier::collapseZeroSet()
    {
        mTimerZeroSetHalfEdge.start();
        EdgeCollapse collapserHalfEdge(mTriangulation, mShell, EdgeCollapse::ZERO_SET, true, 1e-6);
        collapserHalfEdge.collapse();
        mTimerZeroSetHalfEdge.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeTetra(mOutputDirectory + "/zero_set_collapsed_shell(half_edge).vtk", mTriangulation.vertices, mTriangulation.cells);
            jtf::mesh::save_obj((mOutputDirectory + "/zero_set_collapsed_zero_set(half_edge).obj").c_str(), mTriangulation.getZeroSet().triangles,
                                       mTriangulation.getZeroSet().vertices);
            mTriangulation.outputBoundaryFaces(mOutputDirectory + "/zero_set_collapsed(half_edge)");
        }


        mTimerZeroSetGeneral.start();
        EdgeCollapse collapserGeneral(mTriangulation, mShell, EdgeCollapse::ZERO_SET, false, 0.005);
        collapserGeneral.collapse();
        mTimerZeroSetGeneral.end();

        mTriangulation.buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeTetra(mOutputDirectory + "/zero_set_collapsed_shell(general).vtk", mTriangulation.vertices, mTriangulation.cells);
            jtf::mesh::save_obj((mOutputDirectory + "/zero_set_collapsed_zero_set(general).obj").c_str(), mTriangulation.getZeroSet().triangles,
                                       mTriangulation.getZeroSet().vertices);
            mTriangulation.outputBoundaryFaces(mOutputDirectory + "/zero_set_collapsed(general)");
        }
    }

    void Simplifier::writeSummary()
    {
        logTimer("Refinement : ", mTimerRefine);
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
