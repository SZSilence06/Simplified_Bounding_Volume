#ifndef WKY_REFINEMENT_H
#define WKY_REFINEMENT_H

#include "Common.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

namespace SBV
{
    class Refinement
    {
    public:
        Refinement(const matrixr_t& innerShell, const matrixr_t& outerShell);

        bool refine(std::vector<size_t>& output_refinement);

    private:
        enum PointType
        {
            POINT_BOUNDING_BOX,
            POINT_INNER,
            POINT_OUTER
        };

        struct PointInfo
        {
            PointType pointType;
            size_t index;    //index in the original shell samples.
        };

        struct FaceInfo
        {
            std::vector<PointInfo> points;   //points inside the cell.
        };

        using K = CGAL::Exact_predicates_inexact_constructions_kernel;
        using VertexBase = CGAL::Triangulation_vertex_base_with_info_3<PointInfo, K>;
        using CellBase = CGAL::Triangulation_cell_base_with_info_3<FaceInfo, K>;
        using Tds = CGAL::Triangulation_data_structure_3<VertexBase, CellBase, CGAL::Fast_location>;
        using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
        using Point = Delaunay::Point;

    private:
        void init();
        void computeBoundingBox();
        void initErrors();
        bool isFinished();

    private:
        const matrixr_t& mInnerShell;
        const matrixr_t& mOuterShell;

        std::vector<double> mInnerError;    //error for samples on inner shell
        std::vector<double> mOuterError;    //error for samples on outer shell

        Delaunay mDelaunay;
    };
}

#endif
