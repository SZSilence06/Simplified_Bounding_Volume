#include "FMM.h"
#include <cmath>
#include <complex>
#include <boost/math/special_functions/legendre.hpp>
#include <iostream>

using namespace zjucad::matrix;

namespace SBV
{
    int g_factorial[9];
    double g_factorial_reciprocal[9];

    inline void precomputeFactorial()
    {
        g_factorial[0] = 1;
        for(int i = 1; i < 9; i++)
            g_factorial[i] = g_factorial[i - 1] * i;
        for(int i = 1; i < 9; i++)
            g_factorial_reciprocal[i] = 1.0 / g_factorial[1];
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

    inline static std::complex<double> R(const vec3_t& x, int m, int n)
    {
        double r, theta, phi;
        toSphericalCoordinate(x, r, theta, phi);
        double coeff = g_factorial_reciprocal[n + m] * boost::math::legendre_p(n, m, std::cos(theta)) * std::pow(r, n);
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

    void FMM::build(const matrixr_t &vertices, const matrixr_t &triangles, const std::vector<double> &boundary_derivatives)
    {
        precomputeFactorial();
        allocateCells();
        computeFinestStep(vertices);
        initFinestLevel(vertices, triangles, boundary_derivatives);
        upwardPass();
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

    void FMM::computeFinestStep(const matrixr_t &vertices)
    {
        computeBBox(vertices);

        int maxRes = std::pow(2, mMaxLevel);
        mXStep = (mXMax - mXMin) / (maxRes - 1);
        mYStep = (mYMax - mYMin) / (maxRes - 1);
        mZStep = (mZMax - mZMin) / (maxRes - 1);
    }

    void FMM::computeBBox(const matrixr_t &vertices)
    {
        mXMin = std::numeric_limits<double>::max();
        mXMax = std::numeric_limits<double>::lowest();
        mYMin = std::numeric_limits<double>::max();
        mYMax = std::numeric_limits<double>::lowest();
        mZMin = std::numeric_limits<double>::max();
        mZMax = std::numeric_limits<double>::lowest();
        for(size_t i = 0; i < vertices.size(2); i++) {
            double x = vertices(0, i);
            double y = vertices(1, i);
            double z = vertices(2, i);
            mXMin = mXMin < x ? mXMin : x;
            mXMax = mXMax > x ? mXMax : x;
            mYMin = mYMin < y ? mYMin : y;
            mYMax = mYMax > y ? mYMax : y;
            mZMin = mZMin < z ? mZMin : z;
            mZMax = mZMax > z ? mZMax : z;
        }
    }

    void FMM::initFinestLevel(const matrixr_t &vertices, const matrixr_t &triangles, const std::vector<double> &boundary_derivatives)
    {
        for(size_t i = 0; i < triangles.size(2); i++) {
            FMMFace face;
            face.triangle(colon(), 0) = vertices(colon(), triangles(0, i));
            face.triangle(colon(), 1) = vertices(colon(), triangles(1, i));
            face.triangle(colon(), 2) = vertices(colon(), triangles(2, i));
            face.derivative = boundary_derivatives[i];

            vec3_t center = centroid(face.triangle);
            size_t xIndex, yIndex, zIndex;
            getCellIndex(center, xIndex, yIndex, zIndex);
            CellPtr cell = mCells[mMaxLevel - 1][xIndex][yIndex][zIndex];
            if(cell == nullptr) {
                cell = createCell(mMaxLevel - 1, xIndex, yIndex, zIndex);
            }
            cell->faces.push_back(face);
        }
    }

    void FMM::getCellIndex(const vec3_t &center, size_t &xIndex, size_t &yIndex, size_t &zIndex)
    {
        xIndex = (center[0] - mXMin) / (mXStep);
        yIndex = (center[1] - mYMin) / (mYStep);
        zIndex = (center[2] - mZMin) / (mZStep);
    }

    CellPtr FMM::createCell(size_t level, size_t xIndex, size_t yIndex, size_t zIndex)
    {
        CellPtr cell = CellPtr(new FMMCell());
        mCells[level][xIndex][yIndex][zIndex] = cell;
        cell->centroid[0] = mXMin + (xIndex + 0.5) * mXStep;
        cell->centroid[1] = mYMin + (yIndex + 0.5) * mYStep;
        cell->centroid[2] = mZMin + (zIndex + 0.5) * mZStep;
        mActiveCells[level].push_back(cell);

        // handle parent cells
        if(level > 0) {
            size_t xIndexParent = xIndex / 2;
            size_t yIndexParent = yIndex / 2;
            size_t zIndexParent = zIndex / 2;
            CellPtr parentCell = mCells[level - 1][xIndexParent][yIndexParent][zIndexParent];
            if(cell == nullptr) {
                parentCell = createCell(level - 1, xIndexParent, yIndexParent, zIndexParent);
            }
            cell->parent = parentCell;
            parentCell->children.push_back(cell);
        }

        return cell;
    }

    void FMM::upwardPass()
    {
        computeMultipoleForFinest();
    }

    void FMM::computeMultipoleForFinest()
    {
        auto& active_finest_cells = mActiveCells[mMaxLevel - 1];
        for(auto& cell : active_finest_cells) {
            for(int n = 0; n <= mMaxLevel; n++) {
                for(int m = -n; m <= n; m++) {
                    cell->multipoleExp.moment[n][m+n].real(0);
                    cell->multipoleExp.moment[n][m+n].imag(0);
                    for(auto& face: cell->faces) {
                        std::complex<double> result;
                        computeMultipoleMoment(face, cell->centroid, m, n, result);
                        cell->multipoleExp.moment[n][m+n] += result;
                    }
                }
            }
        }
    }

    void FMM::computeMultipoleMoment(const FMMFace &face, const vec3_t &xc, int m, int n, std::complex<double>& result)
    {
        int count = integrateR(face.triangle, xc, m, n, result);
        std::cout << "integrate count : " << count << std::endl;
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
}
