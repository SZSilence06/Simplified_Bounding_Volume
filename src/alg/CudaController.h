#ifndef WKY_CUDA_CONTROLLER_H
#define WKY_CUDA_CONTROLLER_H

#include "Common.h"
#include "eigen3.3/Eigen/Dense"

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
        
        bool findCollapsePos_Boundary(bool isInner, const Eigen::Matrix3d& Q1, const Eigen::Matrix3d& Q2,
                                      matrixr_t& position, double& out_error);
        bool findCollapsePos_ZeroSet(const Eigen::Matrix3d& Q1, const Eigen::Matrix3d& Q2, double sampleRadius,
                                     matrixr_t& position, double& out_error);
    private:
        CudaControllerImpl* impl = nullptr;
    };
}

#endif
