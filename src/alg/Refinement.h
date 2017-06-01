#ifndef WKY_REFINEMENT_H
#define WKY_REFINEMENT_H

#include "Common.h"
#include "TriangulatedShell.h"
#include "Shell.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#ifdef VER_2D
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#else
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#endif

namespace SBV
{
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
            PointInfo(PointType pointType, int index) : pointType(pointType), index(index) {}
            bool operator == (const PointInfo& info) const
            {
                return this->pointType == info.pointType
                        && this->index == info.index;
            }
        };

        struct FaceInfo
        {
            std::vector<PointInfo> points;   //points inside the cell.
            PointInfo maxErrorPoint;
            PointInfo v0, v1, v2;         //vertices of the face, used for judging whether the face is new
        };

#ifdef VER_2D
        using K = CGAL::Exact_predicates_inexact_constructions_kernel;
        using VertexBase = CGAL::Triangulation_vertex_base_with_info_2<PointInfo, K>;
        using CellBase = CGAL::Triangulation_face_base_with_info_2<FaceInfo, K>;
        using Tds = CGAL::Triangulation_data_structure_2<VertexBase, CellBase>;
        using Delaunay = CGAL::Delaunay_triangulation_2<K, Tds>;
        using DPoint = Delaunay::Point;
        using Cell = Delaunay::Face;
        using VertexHandle = Delaunay::Vertex_handle;
#else
        using K = CGAL::Exact_predicates_inexact_constructions_kernel;
        using VertexBase = CGAL::Triangulation_vertex_base_with_info_3<PointInfo, K>;
        using CellBase = CGAL::Triangulation_cell_base_with_info_3<FaceInfo, K>;
        using Tds = CGAL::Triangulation_data_structure_3<VertexBase, CellBase>;
        using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
        using DPoint = Delaunay::Point;
        using Cell = Delaunay::Cell;
        using VertexHandle = Delaunay::Vertex_handle;
#endif

    private:
        void init();
        void computeBoundingBox();
        void initErrors();
        void updateErrors();
        void updatePointInCell(Cell& cell);
        double computeFValue(const Point& point, const Cell& cell);
        double getFValue(const VertexHandle& vh);
        double getError(const PointInfo& point);
        void getPointMatrix(const PointInfo& point, Point& pointMatrix);
        bool isFinished();
        bool isNewCell(const Cell& cell);
        bool checkCondition3(const Cell& cell);
        bool checkClassification(const Cell& cell, const Point& point, bool isOuter);
        double computeHeight(const Cell& cell);
#ifdef VER_2D
        void computeAABB(const Cell& cell, double& xmax, double& xmin, double& ymax, double& ymin);
#else
        void computeAABB(const Cell& cell, double& xmax, double& xmin, double& ymax, double& ymin, double& zmax, double& zmin);
#endif

    private:
        const Shell& mShell;
        TriangulatedShell& mOutput;
        double mAlpha = 0.2;
        double mSampleRadius;

        std::vector<double> mInnerError;    //error for samples on inner shell
        std::vector<double> mOuterError;    //error for samples on outer shell

        PointInfo mNextInsertPoint;         //next point to insert(with maximum error)

        Delaunay mDelaunay;
    };
}

#endif
