#ifndef WKY_CUDA_SHELL_H
#define WKY_CUDA_SHELL_H

#include "eigen3/Eigen/Dense"

namespace SBV
{
    class CudaShell
    {
    public:
        Eigen::MatrixXd innerShell;
        Eigen::MatrixXd outerShell;
    };
}

#endif
