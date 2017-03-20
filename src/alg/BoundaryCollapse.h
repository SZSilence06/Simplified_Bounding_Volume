#ifndef WKY_BOUNDARY_COLLAPSE_H
#define WKY_BOUNDARY_COLLAPSE_H

#include "Refinement.h"
#include <queue>
#include <eigen3/Eigen/Dense>

namespace SBV
{
    class BoundaryCollapse
    {
    public:
        BoundaryCollapse(TriangulatedShell& triangulation, const matrixr_t& innerShell, const matrixr_t& outerShell);

        void collapse();

    private:
        struct EdgeInfo
        {
            size_t firstVert;          //the vert to collapse
            size_t secondVert;         //the vert to collapse to
            double error;         //error for collapsing firstVert to secondVert
            bool isUpdated = false;         // indicating whether this info is out-dated(need to be discarded)
        };

        struct Compare
        {
            bool operator()(std::shared_ptr<EdgeInfo> edge1, std::shared_ptr<EdgeInfo> edge2)
            {
                return edge1->error < edge2->error;
            }
        };

        std::priority_queue<std::shared_ptr<EdgeInfo>, std::vector<std::shared_ptr<EdgeInfo> >, Compare> mQueue;

    private:
        void buildEdgeInfo();
        void buildMatrices();
        void computeErrorMatrix(size_t vert);
        void addNeighbour(size_t firstVert, size_t secondVert);
        void tryAddCollapseableNeighbour(size_t firstVert, size_t secondVert);
        double computeError(size_t vert, const matrixr_t& point);
        bool isValidCollapse(size_t firstVert, size_t secondVert, const matrixr_t &collapseTo);
        void collapseEdge(size_t firstVert, size_t secondVert, const matrixr_t &collapseTo);
        void updateEdgeInfo(size_t vertCollapsed, size_t vertCollapsedTo);
        void insertEdges(size_t vert);
        void mergeNeighbours(size_t vertCollapsed, size_t vertCollapsedTo);
        void organizeOutput();
        size_t getCollapsedVert(size_t vert);
        bool testLinkCondition(size_t firstVert, size_t secondVert);
        void buildOneRingArea(size_t firstVert, size_t secondVert, matrixs_t& lines);
        void findBoundaryEdge(size_t firstVert, size_t secondVert, std::vector<std::pair<size_t, size_t>>& boundaryEdges);

    private:
        TriangulatedShell& mTriangulation;
        const matrixr_t& mInnerShell;
        const matrixr_t& mOuterShell;

        std::vector<Eigen::Matrix3d> mQ;
        std::vector<std::set<size_t> > mNeighbours;   //store the neighbour vertices for each vertex
        std::vector<std::set<size_t> > mCollapseableNeighbours;   //store the collapseable neighbour vertices for each vertex
        std::vector<std::set<size_t> > mNeighbourFaces;

        //for half edge collapse
        std::vector<size_t> mCollapseTo;    //record what vertex the vertices has collapsed to.
        std::vector<std::vector<std::shared_ptr<EdgeInfo> > > mRelatedEdgeInfo;    //recording the edge infos related to the vertices

    };
}

#endif
