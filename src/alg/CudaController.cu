#include "Common.h"
#include "CudaController.h"
#include "CudaControllerImpl.h"
#include "CudaSamplingTree.h"
#include "CudaTriangulatedShell.h"
#include "TriangulatedShell.h"
#include "KernelRegion.h"
#include <wkylib/Geometry/Util.h>

using namespace WKYLIB::Geometry;

namespace SBV
{
    CudaController::~CudaController()
    {
        if(impl)
        {
            delete impl;
        }
    }

    void CudaController::build(const matrixr_t &innerShell, const matrixr_t &outerShell, const TriangulatedShell &triangulation)
    {
        this->impl = new CudaControllerImpl();
        this->impl->build(innerShell, outerShell, triangulation);
    }

    void CudaController::sample(double xmin, double xmax, double ymin, double ymax, double sampleRadius,
                                std::vector<matrixr_t> &output_samples)
    {
        this->impl->sample(xmin, xmax, ymin, ymax, sampleRadius, output_samples);
    }

    void CudaController::buildKernelRegion(const KernelRegion &kernel)
    {
        this->impl->buildKernelRegion(kernel);
    }

    //////////////////////////////////////////////////////////////////////////////////
    void CudaControllerImpl::build(const matrixr_t &innerShell, const matrixr_t &outerShell, const TriangulatedShell &triangulation)
    {
        mShell.assign(CudaShell());
        zju_mat_to_eigen(innerShell, mShell->innerShell);
        zju_mat_to_eigen(outerShell, mShell->outerShell);
        castTriangulation(triangulation, mTriangulation);
    }

    void CudaControllerImpl::castTriangulation(const TriangulatedShell &triangulation, CudaPointer<CudaTriangulatedShell> &cuda_triangulation)
    {
        cuda_triangulation.assign(CudaTriangulatedShell());
        zju_mat_to_eigen(triangulation.vertices, cuda_triangulation->vertices);
        zju_mat_to_eigen(triangulation.triangles, cuda_triangulation->triangles);
        cuda_triangulation->vertType.assign(triangulation.vertType);
    }

    void CudaControllerImpl::buildKernelRegion(const KernelRegion &kernel)
    {
        mKernel.assign(CudaKernelRegion(kernel, mShell, mTriangulation));
    }

    void CudaControllerImpl::sample(double xmin, double xmax, double ymin, double ymax, double sampleRadius,
                                    std::vector<matrixr_t> &output_samples)
    {
        CudaSamplingTree tree(mKernel, xmax, xmin, ymax, ymin, sampleRadius);

        output_samples.reserve(tree.getSamples().size());
        for(int i = 0; i < tree.getSamples().size(); i++)
        {
            matrixr_t mat;
            eigen_to_zju_vector2d(tree.getSamples()[i], mat);
            output_samples.push_back(mat);
        }
    }
}