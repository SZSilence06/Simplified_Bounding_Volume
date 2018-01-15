#include "FMM.h"
#include <cmath>
#include <complex>
#include <boost/math/special_functions/legendre.hpp>
#include <iostream>
#include <cstring>

using namespace zjucad::matrix;

namespace SBV
{
    int g_factorial[9];
    double g_factorial_reciprocal[9];

    inline void precomputeFactorial()
    {
        g_factorial[0] = g_factorial_reciprocal[0] = 1;
        for(int i = 1; i < 9; i++)
            g_factorial[i] = g_factorial[i - 1] * i;
        for(int i = 1; i < 9; i++)
            g_factorial_reciprocal[i] = 1.0 / g_factorial[i];
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

    inline static std::complex<double> S(const vec3_t x, int m, int n)
    {
        double r, theta, phi;
        toSphericalCoordinate(x, r, theta, phi);
        double coeff = g_factorial[n - m] * boost::math::legendre_p(n, m, std::cos(theta)) / std::pow(r, n + 1);
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
        double l1 = sqrt((l3-u3) * (l3-u3) + v3*v3);
        double l2 = sqrt(u3*u3 + v3*v3);

        // threshold for small numbers
        double threshold = 1e-6 * std::min(std::min(l1,l2), l3);
        if(fabs(w0) < threshold)
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
        tplus[0] = sqrt((u3-u0)*(u3-u0) + (v3-v0)*(v3-v0));
        tplus[1] = sqrt(u0*u0 + v0*v0);
        tplus[2] = sqrt((l3-u0)*(l3-u0) + v0*v0);
        tminus[0] = tplus[2];
        tminus[1] = tplus[0];
        tminus[2] = tplus[1];

        // line 1, pp. 1450
        vec3_t R0;
        for(int i = 0; i < 3; i++)
            R0[i] = sqrt(t0[i]*t0[i] + w0*w0);

        //line 2, pp. 1450
        vec3_t Rplus, Rminus;
        for(int i = 0; i < 3; i++)
        {
            Rplus[i] = sqrt(tplus[i]*tplus[i] + w0*w0);
            Rminus[i] = sqrt(tminus[i]*tminus[i] + w0*w0);
        }

        // eq (11)
        vec3_t f2;
        for(int i = 0; i < 3; i++)
        {
            double temp;
            if(w0 == 0)
            {
                if(fabs(t0[i]) < threshold)
                    temp = fabs(log(splus[i]) / sminus[i]);
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
            f2[i] = log(temp);
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
                if(fabs(t0[i]) < threshold)
                    beta[i] = 0;
                else
                    beta[i] = atan(splus[i] / t0[i]) - atan(sminus[i] / t0[i]);
            }
        }
        else
        {
            for(int i = 0; i < 3; i++)
                beta[i] = atan((t0[i]*splus[i]) / (R0[i]*R0[i] + Rplus[i]*fabs(w0))) - atan((t0[i]*sminus[i]) / (R0[i]*R0[i] + Rminus[i]*fabs(w0)));
        }
        betaSum = beta[0] + beta[1] + beta[2];


        // eq (19), integral of kernel 1/R
        double I1 = 0;
        for(int i = 0; i < 3; i++)
            I1 += t0[i]*f2[i];
        I1 -= fabs(w0) * betaSum;
        return I1;
    }

    void FMM::build(const std::vector<mat3x3_t>& triangles, const std::vector<double> &boundary_derivatives)
    {
        precomputeFactorial();
        allocateCells();
        computeFinestStep(triangles);
        initFinestLevel(triangles, boundary_derivatives);
        upwardPass();
        downwardPass();
    }

    double FMM::getPotential(const vec3_t &x)
    {
        size_t xIndex, yIndex, zIndex;
        getCellIndex(x, xIndex, yIndex, zIndex, mMaxLevel - 1);
        double result = 0;

        CellPtr cell = mCells[mMaxLevel - 1][xIndex][yIndex][zIndex];
        if(cell == nullptr) {
            cell = createCell(mMaxLevel - 1, xIndex, yIndex, zIndex, false);
        }

        CellPtr cell2 = cell;
        // evaluate using multipole expansion
        while(cell2 != nullptr) {
            for(auto& interCell : cell2->interList) {
                result += multipoleEvaluate(interCell->multipoleExp, x, interCell->centroid);
            }
            cell2 = cell2->parent;
        }

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

        mXStep.resize(mMaxLevel);
        mYStep.resize(mMaxLevel);
        mZStep.resize(mMaxLevel);

        mXStep[mMaxLevel - 1] = (mXMax - mXMin) / maxRes;
        mYStep[mMaxLevel - 1] = (mYMax - mYMin) / maxRes;
        mZStep[mMaxLevel - 1] = (mZMax - mZMin) / maxRes;
        for(int level = mMaxLevel - 2; level >= 0; level--) {
            mXStep[level] = mXStep[level+1] * 2;
            mYStep[level] = mYStep[level+1] * 2;
            mZStep[level] = mZStep[level+1] * 2;
        }
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
        if(fabs(center[0] - mXMin) < 1e-6)
            xIndex = 0;
        else if(fabs(center[0] - mXMax) < 1e-6)
            xIndex = std::pow(2, level) - 1;
        else
            xIndex = (center[0] - mXMin) / (mXStep[level]);

        if(fabs(center[1] - mYMin) < 1e-6)
            yIndex = 0;
        else if(fabs(center[1] - mYMax) < 1e-6)
            yIndex = std::pow(2, level) - 1;
        else
            yIndex = (center[1] - mYMin) / (mYStep[level]);

        if(fabs(center[2] - mZMin) < 1e-6)
            zIndex = 0;
        else if(fabs(center[2] - mZMax) < 1e-6)
            zIndex = std::pow(2, level) - 1;
        else
            zIndex = (center[2] - mZMin) / (mZStep[level]);
    }

    CellPtr FMM::createCell(size_t level, size_t xIndex, size_t yIndex, size_t zIndex, bool hasFace)
    {
        CellPtr cell = CellPtr(new FMMCell());
        mCells[level][xIndex][yIndex][zIndex] = cell;
        cell->centroid[0] = mXMin + (xIndex + 0.5) * mXStep[level];
        cell->centroid[1] = mYMin + (yIndex + 0.5) * mYStep[level];
        cell->centroid[2] = mZMin + (zIndex + 0.5) * mZStep[level];
        cell->xIndex = xIndex;
        cell->yIndex = yIndex;
        cell->zIndex = zIndex;
        cell->hasFace = hasFace;
        if(hasFace)
            mActiveCells[level].push_back(cell);

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
                    for(int n = 0; n <= mOrder; n++) {
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

    }

    void FMM::computeMultipoleForFinest()
    {
        auto& active_finest_cells = mActiveCells[mMaxLevel - 1];
        for(auto& cell : active_finest_cells) {
            for(int n = 0; n <= mOrder; n++) {
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
        for(int n = 0; n <= mOrder; n++) {
            for(int m = -n; m <= n; m++) {
                for(int n2 = 0; n2 <= n; n2++) {
                    for(int m2 = -n2; m2 <= n2; m2++) {
                        result.moment(n,m) += R(xc - xc_2, m2, n2) * inputMoment.moment(n - n2,m - m2);
                    }
                }
            }
        }
    }

    double FMM::multipoleEvaluate(const MultipoleExp &mul, const vec3_t &x, const vec3_t &xc)
    {
        std::complex<double> sum;
        for(int n = 0; n <= mOrder; n++) {
            for(int m = -n; m <= n; m++) {
                sum += std::conj(S(x - xc, m, n)) * mul.moment(n, m);
            }
        }
        return sum.real();
    }

    inline double FMM::directEvaluate(const FMMFace &face, const vec3_t &x)
    {
        return integrateOverTriangle(x, face.triangle) * face.derivative;
    }
}
