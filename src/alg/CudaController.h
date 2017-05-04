#ifndef WKY_CUDA_CONTROLLER_H
#define WKY_CUDA_CONTROLLER_H

#include "Common.h"

namespace SBV
{
    class CudaControllerImpl;
    class TriangulatedShell;
    class KernelRegion;

    class CudaController
    {
    public:
        CudaController() = default;
        ~CudaController();

        void build(const matrixr_t& innerShell, const matrixr_t& outerShell, const TriangulatedShell& triangulation);

        void sample(double xmin, double xmax, double ymin, double ymax, double sampleRadius,
                    std::vector<matrixr_t>& output_samples);

        void buildKernelRegion(const KernelRegion& kernel);

    private:
        CudaControllerImpl* impl = nullptr;
    };
}

#endif
