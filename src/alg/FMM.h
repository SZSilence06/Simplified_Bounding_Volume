#ifndef WKY_SBV_FMM_H
#define WKY_SBV_FMM_H

#include "Common.h"
#include <eigen3/Eigen/Dense>
#include <complex>
#include <mutex>

namespace SBV
{
    const int MAX_ORDER = 4;

    class FMMCell;
    using CellPtr = std::shared_ptr<FMMCell>;

    class MultipoleExp
    {
    public:
        std::complex<double>& moment(int n, int m) { return _moment[n][m + MAX_ORDER]; }
        const std::complex<double>& moment(int n, int m) const { return _moment[n][m + MAX_ORDER]; }
    private:
        std::complex<double> _moment[MAX_ORDER+1][MAX_ORDER * 2 + 1];
    };

    class LocalExp
    {
    public:
        std::complex<double>& moment(int n, int m) { return _moment[n][m + MAX_ORDER]; }
        const std::complex<double>& moment(int n, int m) const { return _moment[n][m + MAX_ORDER]; }
    private:
        std::complex<double> _moment[MAX_ORDER+1][MAX_ORDER * 2 + 1];
    };

    struct FMMFace
    {
        mat3x3_t triangle;
        double derivative;       // charge density
    };

    struct FMMCell
    {
        vec3_t centroid;                  /// centroid of the cell
        MultipoleExp multipoleExp;        /// multipole expansion of the cell
        LocalExp localExp;                /// local expansion of the cell
        std::vector<FMMFace> faces;       /// faces contained in the cell
        CellPtr parent = nullptr;         /// parent cell
        std::vector<CellPtr> children;    /// children cells
        std::vector<CellPtr> interList;   /// interaction list
        std::vector<CellPtr> neighbours;  /// neighbour cells
        int xIndex, yIndex, zIndex;       /// index of the cell
        int level;                        /// level of the cell
        bool hasFace;                     /// indicating whether there is triangle in this cell
    };

    /**
     * @brief The FMM class
     *
     * This class is used to perform FMM field computation on CPU.
     *
     * Currently, the algorithm seems to have some bugs, causing it to produce low accuracy.
     * One possible reason is that when building the octree, we assign triangles to cells according to the barycenter of the triangle.
     * This may lead to wrong classification for large triangles that cross more than one cell.
     *
     * One possible fix is to slice those triangles that cross cells into small pieces, and make sure each piece is completely inside a cell.
     */
    class FMM
    {
    public:
        FMM() = default;

        double getPotential(const vec3_t& x);

        void setMaxLevel(int maxLevel) { mMaxLevel = maxLevel; }

        void build(const std::vector<mat3x3_t>& triangles, const std::vector<double>& boundary_derivatives);

        //double testGPU(const vec3_t& x);

    private:
        void allocateCells();
        void computeBBox(const std::vector<mat3x3_t> &triangles);
        void computeFinestStep(const std::vector<mat3x3_t>& triangles);
        void initFinestLevel(const std::vector<mat3x3_t>& triangles, const std::vector<double> &boundary_derivatives);
        void getCellIndex(const vec3_t& center, size_t& xIndex, size_t& yIndex, size_t& zIndex, size_t level);
        CellPtr createCell(size_t level, size_t xIndex, size_t yIndex, size_t zIndex, bool hasFace = true);
        void upwardPass();
        void downwardPass();
        void computeMultipoleForFinest();
        void computeMultipoleMoment(const FMMFace &face, const vec3_t &xc, int m, int n, std::complex<double>& result);
        int integrateR(const mat3x3_t& triangle, const vec3_t& xc, int m, int n, std::complex<double>& result);
        void M2M(const MultipoleExp& inputMoment, const vec3_t& xc, const vec3_t& xc_2, MultipoleExp& result);
        void M2L(const MultipoleExp& inputMoment, const vec3_t& xc, const vec3_t& x0, LocalExp& result);
        void L2L(const LocalExp& inputMoment, const vec3_t& x0, const vec3_t& x1, LocalExp& result);
        double multipoleEvaluate(const MultipoleExp& mul, const vec3_t& x, const vec3_t& xc);
        double localEvaluate(const LocalExp& mul, const vec3_t& x, const vec3_t& x0);
        double directEvaluate(const FMMFace& face, const vec3_t& x);

    private:
        int mMaxLevel = 10;

        std::vector<double> mStep;
        double mXMax, mXMin, mYMax, mYMin, mZMax, mZMin;

        std::vector<std::vector<std::vector<std::vector<CellPtr>>>> mCells;
        std::vector<std::vector<CellPtr>> mActiveCells;

        std::recursive_mutex mutex_createCell;

        friend class GPU_FMM;
    };
}

#endif
