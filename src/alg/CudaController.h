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

        void build(const std::vector<Point>& innerShell, const std::vector<Point>& outerShell, const TriangulatedShell& triangulation);

        void buildKernelRegion(const KernelRegion& kernel);
        
        bool findCollapsePos_Boundary(bool isInner, const Eigen::Matrix3d& Q1, const Eigen::Matrix3d& Q2,
                                      Point& position, double& out_error);
        bool findCollapsePos_ZeroSet(const Eigen::Matrix3d& Q1, const Eigen::Matrix3d& Q2, double sampleRadius,
                                     Point& position, double& out_error);
    private:
        CudaControllerImpl* impl = nullptr;
    };
}

#endif
