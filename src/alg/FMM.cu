#include "FMM.h"
#include "GPU_FMM.cuh"
#include <cmath>
#include <complex>
#include <boost/math/special_functions/legendre.hpp>
#include <iostream>
#include <cstring>
#include <sstream>
#include <eigen3/Eigen/Dense>
#include <wkylib/Cuda/CudaVector.h>
#include <wkylib/Cuda/CudaPointer.h>

using namespace zjucad::matrix;
using namespace WKYLIB::Cuda;

namespace SBV {
    namespace __FMM_Internal {  
        static double g_factorial[100];
        static double g_factorial_reciprocal[100];
        static double g_double_factorial[100];

        __global__ void kernel_downwardPass(GPU_FMM* fmm, int* level)
        {
            int x = blockIdx.x * blockDim.x + threadIdx.x;
            int y = blockIdx.y * blockDim.y + threadIdx.y;
            int z = blockIdx.z * blockDim.z + threadIdx.z;

            int res = (int)pow(2.0, (double)*level);
            while(z < res) {
                while(y < res) {
                    while(x < res) {
                        GPU_Cell* cell = fmm->mCells[*level][x][y][z];
                        for(int i = 0; i < cell->interList.size(); i++) {
                            GPU_Cell* interCell = cell->interList[i];
                            GPU_LocalExp localExp;
                            fmm->GPU_M2L(interCell->multipoleExp, interCell->centroid, cell->centroid, localExp);
                            for(int n = 0; n < MAX_ORDER; n++) {
                                for(int m = -n; m <= n; m++) {
                                    cell->localExp.moment(n, m) += localExp.moment(n, m);
                                }
                            }
                        }
                        if(cell->parent) {
                            GPU_LocalExp localExp;
                            fmm->GPU_L2L(cell->parent->localExp, cell->parent->centroid, cell->centroid, localExp);
                            for(int n = 0; n < MAX_ORDER; n++) {
                                for(int m = -n; m <= n; m++) {
                                    cell->localExp.moment(n, m) += localExp.moment(n, m);
                                }
                            }
                        }
                        x += blockDim.x * gridDim.x;
                    }
                    x %= res;
                    y += blockDim.y * gridDim.y;
                }
                y %= res;
                z += blockDim.z * gridDim.z;
            }
        }
    }
}

namespace SBV
{
    using namespace __FMM_Internal;

    inline void precomputeFactorial()
    {
        g_factorial[0] = g_factorial_reciprocal[0] = g_double_factorial[0] = 1;
        for(int i = 1; i < 100; i++)
            g_factorial[i] = g_factorial[i - 1] * i;
        for(int i = 1; i < 100; i++)
            g_factorial_reciprocal[i] = 1.0 / g_factorial[i];
        for(int i = 1; i < 100; i++)
            g_double_factorial[i] = -(2 * i - 1) * g_double_factorial[i - 1];
    }

    inline static vec3_t centroid(const mat3x3_t& triangle)
    {
        return (triangle(colon(), 0) + triangle(colon(), 1) + triangle(colon(), 2)) / 3;
    }

    inline static void toSphericalCoordinate(const vec3_t& x, double& r, double& theta, double& phi)
    {
        r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        theta = std::acos(x[2] / r);
        phi = std::atan2(x[1], x[0]);
    }

    double my_legendre_p(int n, int m, double x)
    {
        if(std::abs(m) > std::abs(n))
            return 0;

        if(m < 0) {
            int sign = (m&1) ? -1 : 1;
            return g_factorial[n+m] / g_factorial[n-m] * sign * my_legendre_p(n, -m, x);
        }

        double p1 = std::pow(1 - x * x, (double)m / 2) * g_double_factorial[m];
        if(m == n)
            return p1;

        int l = m + 1;
        double p2 = (2 * m + 1) * x * p1;
        while(l < n) {
            std::swap(p1, p2);
            p2 = ((2 * l + 1) * x * p1 - (l + m) * p2) / (l + 1 - m);
            l++;
        }
        return p2;
    }

    inline static std::complex<double> R(const vec3_t& x, int m, int n)
    {
        double r, theta, phi;
        toSphericalCoordinate(x, r, theta, phi);
        double coeff;
        if(std::abs(m) > std::abs(n))
            coeff = 0;
        else
            coeff = g_factorial_reciprocal[n + m] * my_legendre_p(n, m, std::cos(theta)) * std::pow(r, n);
        double mphi = m * phi;
        std::complex<double> result;
        result.real(coeff * std::cos(mphi));
        result.imag(coeff * std::sin(mphi));
        return result;
    }


    inline static std::complex<double> S(const vec3_t& x, int m, int n)
    {
        double r, theta, phi;
        toSphericalCoordinate(x, r, theta, phi);
        double coeff = g_factorial[n - m] * my_legendre_p(n, m, std::cos(theta)) / std::pow(r, n + 1);
        double mphi = m * phi;
        std::complex<double> result;
        result.real(coeff * std::cos(mphi));
        result.imag(coeff * std::sin(mphi));
        return result;
    }

    inline static double computeArea(const mat3x3_t& triangle)
    {
        vec3_t a = triangle(colon(), 0);
        vec3_t b = triangle(colon(), 1);
        vec3_t c = triangle(colon(), 2);
        vec3_t ab = b - a;
        vec3_t ac = c - a;
        return 0.5 * (norm(cross(ab,ac)));
    }

    void FMM::build(const std::vector<mat3x3_t>& triangles, const std::vector<double> &boundary_derivatives)
    {
        std::cout << "[INFO] Initializing FMM..." << std::endl;

        precomputeFactorial();
        allocateCells();
        computeFinestStep(triangles);
        initFinestLevel(triangles, boundary_derivatives);

        std::cout << "[INFO] Upward Pass..." << std::endl;
        upwardPass();

        std::cout << "[INFO] Downward Pass..." << std::endl;
        downwardPass();
    }

    double FMM::getPotential(const vec3_t &x)
    {
        size_t xIndex, yIndex, zIndex;
        getCellIndex(x, xIndex, yIndex, zIndex, mMaxLevel - 1);
        double result = 0;

        CellPtr cell = mCells[mMaxLevel - 1][xIndex][yIndex][zIndex];
        if(cell == nullptr) {
            std::cout << "[INFO] creating down cell for evaluation" << std::endl;
            cell = createCell(mMaxLevel - 1, xIndex, yIndex, zIndex, false);
        }

        // evaluate using multipole expansion
        // CellPtr cell2 = cell;
        // while(cell2 != nullptr) {
        //    for(auto& interCell : cell2->interList) {
        //        result += multipoleEvaluate(interCell->multipoleExp, x, interCell->centroid);
        //    }
        //    cell2 = cell2->parent;
        // }

        // evaluate using local expansion
        result += localEvaluate(cell->localExp, x, cell->centroid);

        // evaluate directly
        for(auto& neighbour : cell->neighbours) {
            for(auto& face : neighbour->faces) {
                result += directEvaluate(face, x);
            }
        }

        return result;
    }

    void FMM::allocateCells()
    {
        int levelRes = 1;
        mCells.resize(mMaxLevel);
        mActiveCells.resize(mMaxLevel);
        for(int i = 0; i < mMaxLevel; i++) {
            mCells[i].resize(levelRes);
            for(int j = 0; j < levelRes; j++) {
                mCells[i][j].resize(levelRes);
                for(int k = 0; k < levelRes; k++)
                    mCells[i][j][k].resize(levelRes, nullptr);
            }
            levelRes *= 2;
        }
    }

    void FMM::computeFinestStep(const std::vector<mat3x3_t> &triangles)
    {
        computeBBox(triangles);

        int maxRes = std::pow(2, mMaxLevel - 1);

        mStep.resize(mMaxLevel);

        mStep[mMaxLevel - 1] = (mXMax - mXMin) / maxRes;
        for(int level = mMaxLevel - 2; level >= 0; level--) {
            mStep[level] = mStep[level+1] * 2;
        }
    }

    inline static void scaleInterval(double& min, double& max, double scale)
    {
        double center = (min + max) / 2;
        double half_width = (center - min) * scale;
        min = center - half_width;
        max = center + half_width;
    }

    void FMM::computeBBox(const std::vector<mat3x3_t> &triangles)
    {
        mXMin = std::numeric_limits<double>::max();
        mXMax = std::numeric_limits<double>::lowest();
        mYMin = std::numeric_limits<double>::max();
        mYMax = std::numeric_limits<double>::lowest();
        mZMin = std::numeric_limits<double>::max();
        mZMax = std::numeric_limits<double>::lowest();
        for(size_t i = 0; i < triangles.size(); i++) {
            for(size_t j = 0; j < 3; j++) {
                double x = triangles[i](0, j);
                double y = triangles[i](1, j);
                double z = triangles[i](2, j);
                mXMin = mXMin < x ? mXMin : x;
                mXMax = mXMax > x ? mXMax : x;
                mYMin = mYMin < y ? mYMin : y;
                mYMax = mYMax > y ? mYMax : y;
                mZMin = mZMin < z ? mZMin : z;
                mZMax = mZMax > z ? mZMax : z;
            }
        }

        double maxWidth = std::max(std::max(mXMax - mXMin, mYMax - mYMin), mZMax - mZMin);
        double scaleX = maxWidth / (mXMax - mXMin);
        double scaleY = maxWidth / (mYMax - mYMin);
        double scaleZ = maxWidth / (mZMax - mZMin);
        scaleInterval(mXMin, mXMax, scaleX);
        scaleInterval(mYMin, mYMax, scaleY);
        scaleInterval(mZMin, mZMax, scaleZ);
    }

    void FMM::initFinestLevel(const std::vector<mat3x3_t> &triangles, const std::vector<double> &boundary_derivatives)
    {
        for(size_t i = 0; i < triangles.size(); i++) {
            FMMFace face;
            face.triangle = triangles[i];
            face.derivative = boundary_derivatives[i];

            vec3_t center = centroid(face.triangle);
            size_t xIndex, yIndex, zIndex;
            getCellIndex(center, xIndex, yIndex, zIndex, mMaxLevel - 1);
            CellPtr cell = mCells[mMaxLevel - 1][xIndex][yIndex][zIndex];
            if(cell == nullptr) {
                cell = createCell(mMaxLevel - 1, xIndex, yIndex, zIndex);
            }
            cell->faces.push_back(face);
        }
    }

    void FMM::getCellIndex(const vec3_t &center, size_t &xIndex, size_t &yIndex, size_t &zIndex, size_t level)
    {
        if(std::fabs(center[0] - mXMin) < 1e-6)
            xIndex = 0;
        else if(std::fabs(center[0] - mXMax) < 1e-6)
            xIndex = std::pow(2, level) - 1;
        else
            xIndex = (center[0] - mXMin) / (mStep[level]);

        if(std::fabs(center[1] - mYMin) < 1e-6)
            yIndex = 0;
        else if(std::fabs(center[1] - mYMax) < 1e-6)
            yIndex = std::pow(2, level) - 1;
        else
            yIndex = (center[1] - mYMin) / (mStep[level]);

        if(std::fabs(center[2] - mZMin) < 1e-6)
            zIndex = 0;
        else if(std::fabs(center[2] - mZMax) < 1e-6)
            zIndex = std::pow(2, level) - 1;
        else
            zIndex = (center[2] - mZMin) / (mStep[level]);
    }

    CellPtr FMM::createCell(size_t level, size_t xIndex, size_t yIndex, size_t zIndex, bool hasFace)
    {
        CellPtr cell = CellPtr(new FMMCell());
        mCells[level][xIndex][yIndex][zIndex] = cell;
        cell->centroid[0] = mXMin + (xIndex + 0.5) * mStep[level];
        cell->centroid[1] = mYMin + (yIndex + 0.5) * mStep[level];
        cell->centroid[2] = mZMin + (zIndex + 0.5) * mStep[level];
        cell->xIndex = xIndex;
        cell->yIndex = yIndex;
        cell->zIndex = zIndex;
        cell->level = level;
        cell->hasFace = hasFace;
        if(hasFace)
            mActiveCells[level].push_back(cell);

        this->mutex_createCell.lock();
        // handle parent cells
        if(level > 0) {
            size_t xIndexParent = xIndex / 2;
            size_t yIndexParent = yIndex / 2;
            size_t zIndexParent = zIndex / 2;
            CellPtr parentCell = mCells[level - 1][xIndexParent][yIndexParent][zIndexParent];
            if(parentCell == nullptr) {
                parentCell = createCell(level - 1, xIndexParent, yIndexParent, zIndexParent, hasFace);
            }
            cell->parent = parentCell;
            parentCell->children.push_back(cell);
        }

        // find neighbours
        size_t i = level;
        for(int j = (int)xIndex - 1; j <= (int)xIndex + 1; j++) {
            if(j < 0 || j >= mCells[i].size())
                continue;
            for(int k = (int)yIndex - 1; k <= (int)yIndex + 1; k++) {
                if(k < 0 || k >= mCells[i][j].size())
                    continue;
                for(int l = (int)zIndex - 1; l <= (int)zIndex + 1; l++) {
                    if(l < 0 || l >= mCells[i][j][k].size())
                        continue;
                    if(j == xIndex && k == yIndex && l == zIndex)
                        continue;
                    CellPtr neighbour = mCells[i][j][k][l];
                    if(neighbour != nullptr) {
                        if(neighbour->hasFace)
                            cell->neighbours.push_back(neighbour);
                        if(hasFace)
                            neighbour->neighbours.push_back(cell);
                    }
                }
            }
        }

        // find interaction list
        if(cell->parent != nullptr) {
            CellPtr parentCell = cell->parent;
            for(auto& parentNeighbour : parentCell->neighbours) {
                for(auto& candidate : parentNeighbour->children) {
                    if(std::abs(cell->xIndex - candidate->xIndex) >= 2 || std::abs(cell->yIndex - candidate->yIndex) >= 2
                            || std::abs(cell->zIndex - candidate->zIndex) >= 2) {
                        if(candidate->hasFace)
                            cell->interList.push_back(candidate);
                        if(hasFace)
                            candidate->interList.push_back(cell);
                    }
                }
            }
        }
        this->mutex_createCell.unlock();

        return cell;
    }

    void FMM::upwardPass()
    {
        computeMultipoleForFinest();

        //trace up to sum the multipole moments using M2M
        for(int i = mMaxLevel - 2; i >= 0; i--) {
            auto& active_cells = mActiveCells[i];
            for(auto& cell : active_cells) {
                for(auto& child : cell->children) {
                    MultipoleExp result;
                    M2M(child->multipoleExp, child->centroid, cell->centroid, result);
                    for(int n = 0; n <= MAX_ORDER; n++) {
                        for(int m = -n; m <= n; m++) {
                            cell->multipoleExp.moment(n, m) += result.moment(n, m);
                        }
                    }
                }
            }
        }
    }

    void FMM::downwardPass()
    {
        for(int level = 0; level < mMaxLevel; level++) {
            int res = std::pow(2, level);

#pragma omp parallel for
            for(size_t x = 0; x < res; x++) {
                for(size_t y = 0; y < res; y++) {
                    for(size_t z = 0; z < res; z++) {
                        CellPtr cell = mCells[level][x][y][z];
                        if(cell == nullptr) {
                            cell = createCell(level, x, y, z, false);
                        }
                    }
                }
            }
        }

        GPU_FMM gpu_fmm_tmp;
        CudaPointer<GPU_FMM> gpu_fmm(gpu_fmm_tmp);
        gpu_fmm->buildFromCPU(*this);

        dim3 blocks(4, 4, 4);
        dim3 threads(4, 4, 4);

        for(int level = 0; level < mMaxLevel; level++) {
            std::cout << "[INFO] computing level " << level << std::endl;

            CudaPointer<int> gpu_level(level);
            kernel_downwardPass <<< blocks, threads >>> (gpu_fmm.get(), gpu_level.get());
            cudaDeviceSynchronize();

            int res = std::pow(2, level);

            for(size_t x = 0; x < res; x++) {
                for(size_t y = 0; y < res; y++) {
                    for(size_t z = 0; z < res; z++) {
                        GPU_Cell* gpu_cell = gpu_fmm->mCells[level][x][y][z];
                        CellPtr cell = mCells[level][x][y][z];
                        for(int n = 0; n <= MAX_ORDER; n++) {
                            for(int m = -n; m <= n; m++) {
                                cell->localExp.moment(n, m) = gpu_cell->localExp.moment(n, m);
                            }
                        }
                    }
                }
            }
        } 
    }

    void FMM::computeMultipoleForFinest()
    {
        auto& active_finest_cells = mActiveCells[mMaxLevel - 1];
        for(auto& cell : active_finest_cells) {
            for(int n = 0; n <= MAX_ORDER; n++) {
                for(int m = -n; m <= n; m++) {
                    for(auto& face: cell->faces) {
                        std::complex<double> result;
                        computeMultipoleMoment(face, cell->centroid, m, n, result);
                        cell->multipoleExp.moment(n,m) += result;
                    }
                }
            }
        }
    }

    void FMM::computeMultipoleMoment(const FMMFace &face, const vec3_t &xc, int m, int n, std::complex<double>& result)
    {
        int count = integrateR(face.triangle, xc, m, n, result);
        //std::cout << "integrate count : " << count << std::endl;
        result *= face.derivative;
    }

    int FMM::integrateR(const mat3x3_t &triangle, const vec3_t &xc, int m, int n, std::complex<double>& result)
    {
        double size = computeArea(triangle);
        vec3_t center = centroid(triangle);
        if(size < 1e-2) {
            result = size * R(center - xc, m, n);
            return 1;
        }
        mat3x3_t triangle1 = triangle, triangle2 = triangle, triangle3 = triangle;
        triangle1(colon(), 0) = center;
        triangle2(colon(), 1) = center;
        triangle3(colon(), 2) = center;
        std::complex<double> result1, result2, result3;
        int count1 = integrateR(triangle1, xc, m ,n, result1);
        int count2 = integrateR(triangle2, xc, m ,n, result2);
        int count3 = integrateR(triangle3, xc, m ,n, result3);
        result = result1 + result2 + result3;
        return count1 + count2 + count3;
    }

    void FMM::M2M(const MultipoleExp &inputMoment, const vec3_t &xc, const vec3_t &xc_2, MultipoleExp &result)
    {
        for(int n = 0; n <= MAX_ORDER; n++) {
            for(int m = -n; m <= n; m++) {
                for(int n2 = 0; n2 <= n; n2++) {
                    for(int m2 = -n2; m2 <= n2; m2++) {
                        result.moment(n,m) += R(xc - xc_2, m2, n2) * inputMoment.moment(n - n2,m - m2);
                    }
                }
            }
        }
    }

    void FMM::M2L(const MultipoleExp &inputMoment, const vec3_t &xc, const vec3_t &x0, LocalExp &result)
    {
        for(int n = 0; n <= MAX_ORDER; n++) {
            for(int m = -n; m <= n; m++) {
                for(int n2 = 0; n2 <= MAX_ORDER; n2++) {
                    for(int m2 = -n2; m2 <= n2; m2++) {
                        std::complex<double> temp = std::conj(S(x0 - xc, m + m2, n + n2)) * inputMoment.moment(n2, m2);
                        if(n % 2)
                            temp = -temp;
                        result.moment(n,m) += temp;
                    }
                }
            }
        }
    }

    void FMM::L2L(const LocalExp &inputMoment, const vec3_t &x0, const vec3_t &x1, LocalExp &result)
    {
        for(int n = 0; n <= MAX_ORDER; n++) {
            for(int m = -n; m <= n; m++) {
                for(int n2 = n; n2 <= MAX_ORDER; n2++) {
                    for(int m2 = -n2; m2 <= n2; m2++) {
                        result.moment(n,m) += R(x1 - x0, m2 - m, n2 - n) * inputMoment.moment(n2, m2);
                    }
                }
            }
        }
    }

    double FMM::multipoleEvaluate(const MultipoleExp &mul, const vec3_t &x, const vec3_t &xc)
    {
        std::complex<double> sum;
        for(int n = 0; n <= MAX_ORDER; n++) {
            for(int m = -n; m <= n; m++) {
                sum += std::conj(S(x - xc, m, n)) * mul.moment(n, m);
            }
        }
        return sum.real();
    }

    double FMM::localEvaluate(const LocalExp &mul, const vec3_t &x, const vec3_t &x0)
    {
        std::complex<double> sum;
        for(int n = 0; n <= MAX_ORDER; n++) {
            for(int m = -n; m <= n; m++) {
                sum += R(x - x0, m, n) * mul.moment(n, m);
            }
        }
        return sum.real();
    }

    static void viewTransform(const vec3_t &eye, vec3_t ux, vec3_t uz, mat4x4_t &output)
    {
        uz /= norm(uz);
        ux /= norm(ux);
        const vec3_t uy = cross(uz, ux);

        output(0, 0) = ux[0];
        output(0, 1) = ux[1];
        output(0, 2) = ux[2];
        output(0, 3) = -dot(eye, ux);
        output(1, 0) = uy[0];
        output(1, 1) = uy[1];
        output(1, 2) = uy[2];
        output(1, 3) = -dot(eye, uy);
        output(2, 0) = uz[0];
        output(2, 1) = uz[1];
        output(2, 2) = uz[2];
        output(2, 3) = -dot(eye, uz);
        output(3, 0) = 0;
        output(3, 1) = 0;
        output(3, 2) = 0;
        output(3, 3) = 1;
    }

    static void localTransform(const vec3_t &a, const vec3_t &b, const vec3_t &c, mat4x4_t &output)
    {
        vec3_t ux = b - a;
        vec3_t uz = cross(ux, c - a);
        output.resize(4, 4);
        viewTransform(a, ux, uz, output);
    }

    //closed-form calculation according to Graglia 1993.
    static double integrateOverTriangle(const vec3_t& x, const mat3x3_t &triangle)
    {
        mat4x4_t triangle_4;
        triangle_4(colon(0, 2), colon(0, 2)) = triangle;
        triangle_4(3, colon(0, 2)) = ones<double>(1, 3);
        mat4x4_t transform;
        localTransform(triangle(colon(), 0), triangle(colon(), 1), triangle(colon(), 2), transform);
        mat4x4_t localTriangle = transform * triangle_4;

        double l3 = localTriangle(0, 1);
        double u3 = localTriangle(0, 2);
        double v3 = localTriangle(1, 2);

        if(l3 < 0)
            throw std::runtime_error("l3 < 0.");

        vec4_t tempX;
        tempX(colon(0, 2), colon()) = x;
        tempX[3] = 1;
        vec4_t localX = transform * tempX;
        double u0 = localX[0];
        double v0 = localX[1];
        double w0 = localX[2];

        // edge lengths
        double l1 = std::sqrt((l3-u3) * (l3-u3) + v3*v3);
        double l2 = std::sqrt(u3*u3 + v3*v3);

        // threshold for small numbers
        double threshold = 1e-6 * std::min(std::min(l1,l2), l3);
        if(std::fabs(w0) < threshold)
            w0 = 0;

        // eq (3)
        vec3_t sminus, splus;
        sminus[0] = -((l3-u3)*(l3-u0)+v3*v0)/l1;
        sminus[1] = -(u3*(u3-u0)+v3*(v3-v0))/l2;
        sminus[2] = -u0;
        splus[0] = ((u3-l3)*(u3-u0)+v3*(v3-v0))/l1;
        splus[1] = (u3*u0+v3*v0)/l2;
        splus[2] = l3-u0;

        // eq (4)
        vec3_t t0;
        t0[0] = ((u3-l3)*v0+v3*(l3-u0))/l1;
        t0[1] = (v3*u0-u3*v0)/l2;
        t0[2] = v0;

        // eq (5)
        vec3_t tplus, tminus;
        tplus[0] = std::sqrt((u3-u0)*(u3-u0) + (v3-v0)*(v3-v0));
        tplus[1] = std::sqrt(u0*u0 + v0*v0);
        tplus[2] = std::sqrt((l3-u0)*(l3-u0) + v0*v0);
        tminus[0] = tplus[2];
        tminus[1] = tplus[0];
        tminus[2] = tplus[1];

        // line 1, pp. 1450
        vec3_t R0;
        for(int i = 0; i < 3; i++)
            R0[i] = std::sqrt(t0[i]*t0[i] + w0*w0);

        //line 2, pp. 1450
        vec3_t Rplus, Rminus;
        for(int i = 0; i < 3; i++)
        {
            Rplus[i] = std::sqrt(tplus[i]*tplus[i] + w0*w0);
            Rminus[i] = std::sqrt(tminus[i]*tminus[i] + w0*w0);
        }

        // eq (11)
        vec3_t f2;
        for(int i = 0; i < 3; i++)
        {
            double temp;
            if(w0 == 0)
            {
                if(std::fabs(t0[i]) < threshold)
                    temp = std::fabs(std::log(splus[i]) / sminus[i]);
                else
                    temp = (tplus[i]+splus[i]) / (tminus[i]+sminus[i]);
                if(temp < 0)
                    std::cerr << "[WARNING] computing log of negative number. i = " << i
                              << " tplus[0] = " << tplus[0]
                              << " tminus[0] = " << tminus[0]
                              << " splus[0] = " << splus[0]
                              << " sminus[0] = " << sminus[0]
                              << ". line " << __LINE__ << std::endl;
            }
            else
            {
                 temp = (Rplus[i]+splus[i]) / (Rminus[i]+sminus[i]);
                 if(temp < 0)
                     std::cerr << "[WARNING] computing log of negative number. i = " << i << ". line " << __LINE__ << std::endl;
            }
            f2[i] = std::log(temp);
            //fix value for points on the triangle corners
            if(f2[i] != f2[i])  //nan
                f2[i] = 0;
        }


        // eq (13) and eq (14)
        vec3_t beta;
        double betaSum;
        if(w0 == 0)
        {
            for(int i = 0; i < 3; i++)
            {
                if(std::fabs(t0[i]) < threshold)
                    beta[i] = 0;
                else
                    beta[i] = atan(splus[i] / t0[i]) - atan(sminus[i] / t0[i]);
            }
        }
        else
        {
            for(int i = 0; i < 3; i++)
                beta[i] = atan((t0[i]*splus[i]) / (R0[i]*R0[i] + Rplus[i]*std::fabs(w0))) - atan((t0[i]*sminus[i]) / (R0[i]*R0[i] + Rminus[i]*std::fabs(w0)));
        }
        betaSum = beta[0] + beta[1] + beta[2];


        // eq (19), integral of kernel 1/R
        double I1 = 0;
        for(int i = 0; i < 3; i++)
            I1 += t0[i]*f2[i];
        I1 -= std::fabs(w0) * betaSum;
        return I1;
    }

    inline double FMM::directEvaluate(const FMMFace &face, const vec3_t &x)
    {
        return integrateOverTriangle(x, face.triangle) * face.derivative;
    }

    /*double FMM::testGPU(const vec3_t &x)
    {
        GPU_FMM gpu_fmm_tmp;
        CudaPointer<GPU_FMM> gpu_fmm(gpu_fmm_tmp);
        Eigen::Vector3d xx;
        xx[0] = x[0]; xx[1] = x[1]; xx[2] = x[2];
        gpu_fmm->buildFromCPU(*this);
        return gpu_fmm->getPotential(xx);
    }*/
}
