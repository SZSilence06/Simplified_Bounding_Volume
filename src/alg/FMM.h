#ifndef WKY_SBV_FMM_H
#define WKY_SBV_FMM_H

#include "Common.h"
#include <eigen3/Eigen/Dense>
#include <complex>

namespace SBV
{
    const int MAX_ORDER = 4;

    class FMMCell;
    using CellPtr = std::shared_ptr<FMMCell>;

    struct MultipoleExp
    {
        std::complex<double> moment[MAX_ORDER][MAX_ORDER * 2 + 1];
    };

    struct FMMFace
    {
        mat3x3_t triangle;
        double derivative;
    };

    struct FMMCell
    {
        vec3_t centroid;
        MultipoleExp multipoleExp;
        std::vector<FMMFace> faces;
        CellPtr parent = nullptr;
        std::vector<CellPtr> children;
        std::vector<CellPtr> interList;
    };

    class FMM
    {
    public:
        FMM() = default;

        double getPotential(const vec3_t& x);

        void setMaxLevel(int maxLevel) { mMaxLevel = maxLevel; }
        void setOrder(int order) { mOrder = order; }

        void build(const matrixr_t& vertices, const matrixr_t& triangles, const std::vector<double>& boundary_derivatives);

    private:
        void allocateCells();
        void computeBBox(const matrixr_t &vertices);
        void computeFinestStep(const matrixr_t& vertices);
        void initFinestLevel(const matrixr_t &vertices, const matrixr_t &triangles, const std::vector<double> &boundary_derivatives);
        void getCellIndex(const vec3_t& center, size_t& xIndex, size_t& yIndex, size_t& zIndex);
        CellPtr createCell(size_t level, size_t xIndex, size_t yIndex, size_t zIndex);
        void upwardPass();
        void downwardPass();
        void computeMultipoleForFinest();
        void computeMultipoleMoment(const FMMFace &face, const vec3_t &xc, int m, int n, std::complex<double>& result);
        int integrateR(const mat3x3_t& triangle, const vec3_t& xc, int m, int n, std::complex<double>& result);

    private:
        int mMaxLevel = 10;
        int mOrder = 4;

        double mXStep;
        double mYStep;
        double mZStep;
        double mXMax, mXMin, mYMax, mYMin, mZMax, mZMin;

        std::vector<std::vector<std::vector<std::vector<CellPtr>>>> mCells;
        std::vector<std::vector<CellPtr>> mActiveCells;
    };


}

#endif
