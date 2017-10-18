#ifndef WKY_REFINEMENT_H
#define WKY_REFINEMENT_H

#include "Common.h"
#include "TriangulatedShell.h"
#include "Shell.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

namespace SBV
{
    class BaryComputer;

    class Refinement
    {      
    public:
        Refinement(const Shell& shell, TriangulatedShell &output, double alpha, double sampleRadius);

        bool refine();

    private:
        struct PointInfo
        {
            PointType pointType = POINT_UNKNOWN;
            int index = -1;    //index in the original shell samples.
            int indexInDelaunay = -1;    //index in the delaunay triangulation

            PointInfo() {}
            PointInfo(PointType pointType, size_t index) : pointType(pointType), index(index) {}
            bool operator == (const PointInfo& info) const
            {
                return this->pointType == info.pointType
                        && this->index == info.index;
            }
        };

        struct CellInfo
        {
            std::vector<PointInfo> points;   //points inside the cell.
            PointInfo maxErrorPoint;
            bool isNew = true;
            bool isJudged = false;
            bool isBad = false;
            bool notResolvable = false;
        };

        using K = CGAL::Exact_predicates_inexact_constructions_kernel;
        using VertexBase = CGAL::Triangulation_vertex_base_with_info_3<PointInfo, K>;
        using CellBase = CGAL::Triangulation_cell_base_with_info_3<CellInfo, K>;
        using Tds = CGAL::Triangulation_data_structure_3<VertexBase, CellBase>;
        using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
        using Point = Delaunay::Point;
        using Cell = Delaunay::Cell;
        using VertexHandle = Delaunay::Vertex_handle;

    private:
        void init();
        void computeBoundingBox();
        void initErrors();
        void updateErrors();
        void updatePointInCell(Cell& cell);
        void organizeOutput();
        PointType getPointType(const VertexHandle& vh);
        double getFValue(const VertexHandle& vh);
        double getError(const PointInfo& point);
        void getPointMatrix(const PointInfo& point, matrixr_t& pointMatrix);
        bool isFinished();
        bool isBadCell(Cell& cell);
        bool checkCondition3(Cell& cell);
        bool checkClassification(const Cell& cell, const BaryComputer& baryComputer, const matrixr_t& point, bool isOuter);
        bool resolveBadCell();
        double computeHeight(const Cell& cell);
        void computeAABB(const Cell& cell, double& xmax, double& xmin, double& ymax, double& ymin, double& zmax, double& zmin);
        bool insertPoint(const PointInfo& info, VertexHandle& vertexHandle);

    private:
        const Shell& mShell;
        TriangulatedShell& mOutput;
        double mAlpha = 0.2;
        double mSampleRadius;

        std::vector<double> mInnerError;    //error for samples on inner shell
        std::vector<double> mOuterError;    //error for samples on outer shell
        std::vector<bool> mInnerExists;
        std::vector<bool> mOuterExists;

        PointInfo mNextInsertPoint;         //next point to insert(with maximum error)

        Cell* mBadCell = nullptr;     //the bad cell to be resolved next

        Delaunay mDelaunay;
    };
}

#endif
