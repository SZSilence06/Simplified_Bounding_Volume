#include "GPU_FMM.h"
#include "FMM.h"
#include <thrust/swap.h>

namespace SBV
{
    template<class T1, class T2>
    static void toEigen(const T1& zjumat, T2& out)
    {
        for(size_t i = 0; i < zjumat.size(1); i++)
            for(size_t j = 0; j < zjumat.size(2); j++)
                out(i, j) = zjumat(i, j);
    }

    void GPU_FMM::buildFromCPU(const FMM &fmm)
    {
        this->mMaxLevel = fmm.mMaxLevel;
        this->mDownLevel = fmm.mDownLevel;
        this->mOrder = fmm.mOrder;
        this->mXMax = fmm.mXMax;
        this->mXMin = fmm.mXMin;
        this->mYMax = fmm.mYMax;
        this->mYMin = fmm.mYMin;
        this->mZMax = fmm.mZMax;
        this->mZMin = fmm.mZMin;
        this->mStep.assign(fmm.mStep);

        this->factorial[0] = this->factorial_reciprocal[0] = this->double_factorial[0] = 1;
        for(int i = 1; i < 100; i++)
            this->factorial[i] = this->factorial[i - 1] * i;
        for(int i = 1; i < 100; i++)
            this->factorial_reciprocal[i] = 1.0 / this->factorial[i];
        for(int i = 1; i < 100; i++)
            this->double_factorial[i] = -(2 * i - 1) * this->double_factorial[i - 1];

        cudaMallocManaged(&mCells, sizeof(GPU_Cell****) * mMaxLevel);

        int levelRes = 1;
        for(size_t i = 0; i < mMaxLevel; i++) {
            cudaMallocManaged(&mCells[i], sizeof(GPU_Cell***) * levelRes);
            for(size_t j = 0; j < levelRes; j++) {
                cudaMallocManaged(&mCells[i][j], sizeof(GPU_Cell**) * levelRes);
                for(size_t k = 0; k < levelRes; k++) {
                    cudaMallocManaged(&mCells[i][j][k], sizeof(GPU_Cell*) * levelRes);
                    for(size_t l = 0; l < levelRes; l++) {
                        cudaMallocManaged(&mCells[i][j][k][l], sizeof(GPU_Cell));
                    }
                }
            }
            levelRes *= 2;
        }

        levelRes = 1;
        for(size_t i = 0; i < mMaxLevel; i++) {
            for(size_t j = 0; j < levelRes; j++) {
                for(size_t k = 0; k < levelRes; k++) {
                    for(size_t l = 0; l < levelRes; l++) {
                        cudaMallocManaged(&mCells[i][j][k][l], sizeof(GPU_Cell));

                    }
                }
            }
            levelRes *= 2;
        }

        levelRes = 1;
        for(size_t i = 0; i < mMaxLevel; i++) {
            for(size_t j = 0; j < levelRes; j++) {
                for(size_t k = 0; k < levelRes; k++) {
                    for(size_t l = 0; l < levelRes; l++) {
                        CpuCell2Gpu(fmm.mCells[i][j][k][l], mCells[i][j][k][l]);
                    }
                }
            }
            levelRes *= 2;
        }
    }

    void GPU_FMM::CpuCell2Gpu(CellPtr cpu_cell, GPU_Cell* gpu_cell)
    {
        for(size_t i = 0; i < 3; i++)
            gpu_cell->centroid[i] = cpu_cell->centroid[i];

        CpuCellVec2Gpu(cpu_cell->children, gpu_cell->children);
        CpuCellVec2Gpu(cpu_cell->interList, gpu_cell->interList);
        CpuCellVec2Gpu(cpu_cell->neighbours, gpu_cell->neighbours);

        gpu_cell->level = cpu_cell->level;
        gpu_cell->xIndex = cpu_cell->xIndex;
        gpu_cell->yIndex = cpu_cell->yIndex;
        gpu_cell->zIndex = cpu_cell->zIndex;
        gpu_cell->hasFace = cpu_cell->hasFace;

        // copy faces
        std::vector<GPU_Face> gpu_face_vec;
        for(FMMFace& face : cpu_cell->faces) {
            GPU_Face gpu_face;
            toEigen(face.triangle, gpu_face.triangle);
            gpu_face.derivative = face.derivative;
            gpu_face_vec.push_back(gpu_face);
        }
        gpu_cell->faces.assign(gpu_face_vec);

        if(cpu_cell->parent)
            gpu_cell->parent = mCells[cpu_cell->parent->level][cpu_cell->parent->xIndex][cpu_cell->parent->yIndex][cpu_cell->parent->zIndex];

        for(int n = 0; n <= mOrder; n++) {
            for(int m = -n; m <= n; m++) {
                gpu_cell->multipoleExp.moment(n, m) = cpu_cell->multipoleExp.moment(n, m);
                gpu_cell->localExp.moment(n, m) = cpu_cell->localExp.moment(n, m);
            }
        }
    }

    void GPU_FMM::CpuCellVec2Gpu(const std::vector<CellPtr>& cpu_vec, WKYLIB::Cuda::CudaVector<GPU_Cell*> &gpu_vec)
    {
        std::vector<GPU_Cell*> cpu_tmp_vec;
        for(size_t i = 0; i < cpu_vec.size(); i++) {
            CellPtr child = cpu_vec[i];
            cpu_tmp_vec.push_back(mCells[child->level][child->xIndex][child->yIndex][child->zIndex]);
        }
        gpu_vec.assign(cpu_tmp_vec);
    }

    __device__ double GPU_FMM::GPU_legendre_p(int n, int m, double x)
    {
        if(fabs((double)m) > fabs((double)n))
            return 0;

        double coeff = 1.0;
        if(m < 0) {
            int sign = (m&1) ? -1 : 1;
            coeff = this->factorial[n+m] / this->factorial[n-m] * sign;
            m = -m;
        }

        double p1 = pow(1 - x * x, (double)m / 2) * this->double_factorial[m];
        if(m == n)
            return p1 * coeff;

        int l = m + 1;
        double p2 = (2 * m + 1) * x * p1;
        while(l < n) {
            thrust::swap(p1, p2);
            p2 = ((2 * l + 1) * x * p1 - (l + m) * p2) / (l + 1 - m);
            l++;
        }
        return p2 * coeff;
    }

    __device__ inline void GPU_FMM::GPU_toSphericalCoordinate(const Eigen::Vector3d& x, double& r, double& theta, double& phi)
    {
        r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        theta = acos(x[2] / r);
        phi = atan2(x[1], x[0]);
    }

    __device__ inline thrust::complex<double> GPU_FMM::GPU_R(const Eigen::Vector3d& x, int m, int n)
    {
        double r, theta, phi;
        GPU_toSphericalCoordinate(x, r, theta, phi);
        double coeff;
        if(fabs((double)m) > fabs((double)n))
            coeff = 0;
        else
            coeff = this->factorial_reciprocal[n + m] * GPU_legendre_p(n, m, cos(theta)) * pow(r, n);
        double mphi = m * phi;
        thrust::complex<double> result;
        result.real(coeff * cos(mphi));
        result.imag(coeff * sin(mphi));
        return result;
    }

    __device__ inline thrust::complex<double> GPU_FMM::GPU_S(const Eigen::Vector3d& x, int m, int n)
    {
        double r, theta, phi;
        GPU_toSphericalCoordinate(x, r, theta, phi);
        double coeff = this->factorial[n - m] * GPU_legendre_p(n, m, cos(theta)) / pow(r, n + 1);
        double mphi = m * phi;
        thrust::complex<double> result;
        result.real(coeff * cos(mphi));
        result.imag(coeff * sin(mphi));
        return result;
    }

    __device__ void GPU_FMM::GPU_M2L(const GPU_MultipoleExp &inputMoment, const Eigen::Vector3d &xc, const Eigen::Vector3d &x0, GPU_LocalExp &result)
    {
        for(int n = 0; n <= mOrder; n++) {
            for(int m = -n; m <= n; m++) {
                for(int n2 = 0; n2 <= mOrder; n2++) {
                    for(int m2 = -n2; m2 <= n2; m2++) {
                        thrust::complex<double> temp = thrust::conj(GPU_S(x0 - xc, m + m2, n + n2)) * inputMoment.moment(n2, m2);
                        if(n % 2)
                            temp = -temp;
                        result.moment(n,m) += temp;
                    }
                }
            }
        }
    }

    __device__ void GPU_FMM::GPU_L2L(const GPU_LocalExp &inputMoment, const Eigen::Vector3d &x0, const Eigen::Vector3d &x1, GPU_LocalExp &result)
    {
        for(int n = 0; n <= mOrder; n++) {
            for(int m = -n; m <= n; m++) {
                for(int n2 = n; n2 <= mOrder; n2++) {
                    for(int m2 = -n2; m2 <= n2; m2++) {
                        result.moment(n,m) += GPU_R(x1 - x0, m2 - m, n2 - n) * inputMoment.moment(n2, m2);
                    }
                }
            }
        }
    }

    __host__ __device__ double GPU_FMM::getPotential(const Eigen::Vector3d& x)
    {
        /*size_t xIndex, yIndex, zIndex;
        getCellIndex(x, xIndex, yIndex, zIndex, mDownLevel - 1);
        double result = 0;

        GPU_Cell* cell = mCells[mMaxLevel - 1][xIndex][yIndex][zIndex];

        // evaluate using multipole expansion
        GPU_Cell* cell2 = cell;
        while(cell2 != nullptr) {
            for(auto& interCell : cell2->interList) {
                result += multipoleEvaluate(interCell->multipoleExp, x, interCell->centroid);
            }
            cell2 = cell2->parent;
        }

        // evaluate using local expansion
        //result += localEvaluate(cell->localExp, x, cell->centroid);

        // evaluate directly
        for(auto& neighbour : cell->neighbours) {
            for(auto& face : neighbour->faces) {
                result += directEvaluate(face, x);
            }
        }

        return result;*/
    }
}
