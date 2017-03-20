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

        buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/refined_shell.obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/refined_zero_set.obj", mZeroSet.vertices, mZeroSet.lines);
        }
    }

    void Simplifier::collapseBoundary()
    {
        BoundaryCollapse collapser(mTriangulation, mInnerShell, mOuterShell);
        collapser.collapse();

        buildZeroSet();
        if(mNeedGenTempResult)
        {
            WKYLIB::Mesh::writeMesh2D(mOutputDirectory + "/boundary_collapsed_shell.obj", mTriangulation.vertices, mTriangulation.triangles);
            WKYLIB::Mesh::writeCurve2D(mOutputDirectory + "/boundary_collapsed_zero_set.obj", mZeroSet.vertices, mZeroSet.lines);
        }
    }

    void Simplifier::buildZeroSet()
    {
        std::vector<std::pair<size_t, size_t> > existing_verts;

        //find out the number of faces in the zero set
        int zeroFaceCount = 0;
        for(int i = 0; i < mTriangulation.triangles.size(2); i++)
        {
            size_t v0 = mTriangulation.triangles(0, i);
            size_t v1 = mTriangulation.triangles(1, i);
            size_t v2 = mTriangulation.triangles(2, i);

            if(mTriangulation.getSign(v0) == mTriangulation.getSign(v1) &&
                    mTriangulation.getSign(v0) == mTriangulation.getSign(v2))
            {
                //F value sign of the vertices are same, so no zero-set in this triangle
                continue;
            }

            zeroFaceCount++;
        }

        //build zero face connection
        mZeroSet.lines.resize(2, zeroFaceCount);

        for(int i = 0, j = 0; i < mTriangulation.triangles.size(2); i++)
        {
            size_t v0 = mTriangulation.triangles(0, i);
            size_t v1 = mTriangulation.triangles(1, i);
            size_t v2 = mTriangulation.triangles(2, i);

            if(mTriangulation.getSign(v0) == mTriangulation.getSign(v1) &&
                    mTriangulation.getSign(v0) == mTriangulation.getSign(v2))
            {
                //F value sign of the vertices are same, so no zero-set in this triangle
                continue;
            }

            if(mTriangulation.getSign(v0) == mTriangulation.getSign(v1))
            {
                mZeroSet.lines(0, j) = getZeroPointIndex(v0, v2, existing_verts);
                mZeroSet.lines(1, j) = getZeroPointIndex(v1, v2, existing_verts);
            }
            else if(mTriangulation.getSign(v0) == mTriangulation.getSign(v2))
            {
                mZeroSet.lines(0, j) = getZeroPointIndex(v0, v1, existing_verts);
                mZeroSet.lines(1, j) = getZeroPointIndex(v1, v2, existing_verts);
            }
            else if(mTriangulation.getSign(v1) == mTriangulation.getSign(v2))
            {
                mZeroSet.lines(0, j) = getZeroPointIndex(v0, v1, existing_verts);
                mZeroSet.lines(1, j) = getZeroPointIndex(v0, v2, existing_verts);
            }
            else
            {
                throw std::runtime_error("cannot find valid zero point when building zero sets.");
            }
            j++;
        }

        //build zero point positions
        mZeroSet.vertices.resize(2, existing_verts.size());
        for(int i = 0; i< existing_verts.size(); i++)
        {
            auto& vertPair = existing_verts[i];
            mZeroSet.vertices(colon(), i) = (mTriangulation.vertices(colon(), vertPair.first) +
                    mTriangulation.vertices(colon(), vertPair.second)) / 2;
        }
    }

    size_t Simplifier::getZeroPointIndex(size_t firstVertex, size_t secondVertex,
                                         std::vector<std::pair<size_t, size_t> > &existingVertPairs)
    {
        for(int i = 0; i < existingVertPairs.size(); i++)
        {
            auto& vertPair = existingVertPairs[i];
            if((vertPair.first == firstVertex && vertPair.second == secondVertex)
                    || (vertPair.first == secondVertex && vertPair.second == firstVertex))
            {
                return i;
            }
        }

        existingVertPairs.push_back(std::make_pair(firstVertex, secondVertex));
        return existingVertPairs.size() - 1;
    }
}
