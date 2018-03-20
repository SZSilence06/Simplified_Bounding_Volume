#ifndef WKY_EDGE_COLLAPSE_H
#define WKY_EDGE_COLLAPSE_H

#include <queue>
#include <mutex>
#include <eigen3/Eigen/Dense>
#include "TriangulatedShell.h"
#include "Shell.h"

namespace SBV
{
    /**
    * @brief The EdgeCollapse class
    *
    * This class is used to perform the edge collapse operation for the boundary and zero set faces.
    *
    * So far, the algorithm does not have a half-edge mesh structure, so the collapse operation is virtually performed as the following:
    * we maintain a list recording the vertex that a vertex has been collapsed to for all vertices.
    * For example, if vertex 1 has been collapsed to vertex 2, we record collapseTo[1] = 2.
    * For non half-edge collapse, we modify the position of the collapsed-to vertex.
    * After all edge collapse operations are finished, we reconstruct the vertex and triangle matrices from the above information.
    */
    class EdgeCollapse
    {
    public:
        enum Type{
            BOUNDARY,
            ZERO_SET
        };

    public:
        /**
         * @brief EdgeCollapse constructor
         * @param triangulation : the input delaunay triangulation
         * @param shell : the input sample points
         * @param type : collapse type
         * @param isHalfEdge : true if use half edge collapse, false otherwise
         * @param sampleRadius : used for finding candidate collapse targets for general collapse on zero set.
         *                       This parameter determines the sample radius of sampling the kernel region.
         */
        EdgeCollapse(TriangulatedShell& triangulation, const Shell& shell, Type type, bool isHalfEdge, double sampleRadius);

        /**
         * @brief perform the collapse.
         */
        void collapse();

    private:
        struct EdgeInfo
        {
            size_t firstVert;          //the vert to collapse
            size_t secondVert;         //the vert to collapse to
            vec3_t position;        //the position to collapse to
            double error;         //error for collapsing the edge
            bool isUpdated = false;         // indicating whether this info is out-dated(need to be discarded)
        };

        struct Compare
        {
            bool operator()(std::shared_ptr<EdgeInfo> edge1, std::shared_ptr<EdgeInfo> edge2)
            {
                return edge1->error > edge2->error;
            }
        };

        std::priority_queue<std::shared_ptr<EdgeInfo>, std::vector<std::shared_ptr<EdgeInfo> >, Compare> mQueue;

    private:
        /**
         * @brief build initial priority queue for the edge collapse.
         */
        void buildEdgeInfo();

        /**
         * @brief build the Q matrices for the vertices.
         */
        void buildMatrices();

        /**
         * @brief compute the Q matrix for the given vertex.
         * @param vert : The id of the vertex to compute the Q matrix.
         */
        void computeErrorMatrix(size_t vert);

        /**
         * @brief Adds the neighbour vertices of the secondVert to the firstVert.
         */
        void addNeighbour(size_t firstVert, size_t secondVert);

        /**
         * @brief  Adds the collapseable neighbour vertices of the secondVert to the firstVert.
         *
         * The collapseable neighbour vertices are defined as vertices that share the same point type with the given vertex.
         * i.e, inner shell vertices and inner shell vertices are collapseable neighbour vertices,
         * outer shell vertices and outer shell vertices are collapseable neighbour vertices.
         * This is because an edge collapse only occurs on edges whose vertices share the same point type.
         */
        void tryAddCollapseableNeighbour(size_t firstVert, size_t secondVert);

        /**
         * @brief Adds the collapseable neighbour triangles of the three input vertices.
         *
         * The collapseable neighbour triangles are defined as triangles that are on same side with the given vertex.
         * i.e, for inner shell vertices, collapseable neighbour triangles are triangles that are on inner shell.
         */
        void tryAddCollapseableNeighbourFace(size_t firstVert, size_t secondVert, size_t thirdVert);\

        /**
         * @brief compute the quatric error given the vertex and evaluated point.
         * @param vert : the id of the given vertex.
         * @param point : the point to be evaluated.
         * @return the quadric error.
         */
        double computeError(size_t vert, const matrixr_t& point);

        /**
         * @brief determine whether the given three vertices are on same shell.
         */
        bool isSameShellFace(size_t firstVert, size_t secondVert, size_t thirdVert);

        /**
         * @brief determine whether the two given vertices are neighbour vertices according to the definition defined in tryAddCollapseableNeighbour().
         */
        bool isCollapseable(size_t firstVert, size_t secondVert);

        /**
         * @brief determine whether a candidate collapse is valid.
         *        The criterion contains link condition, classification condition and normal condition.
         *
         * @param firstVert : the first vertex of the collapsed edge
         * @param secondVert : the second vertex of the collapsed edge
         * @param collapseTo : the position that the edge is collapsed to
         * @return true if valid, and false otherwise
         */
        bool isValidCollapse(size_t firstVert, size_t secondVert, const vec3_t &collapseTo);

        /**
         * @brief perform an edge collapse operation.
         *
         *        When doing the operation, we do not immediately modify the position and connection of the triangulation.
         *        We just mark mCollapseTo[first] = secondVert, and update the position of the second vertex.
         *        The actual modification of the triangulation is done after all edge collapses are finished.
         *
         * @param firstVert : the first vertex of the collapsed edge
         * @param secondVert : the second vertex of the collapsed edge
         * @param collapseTo : the position that the edge is collapsed to
         */
        void collapseEdge(size_t firstVert, size_t secondVert, const matrixr_t &collapseTo);

        /**
         * @brief update the neighbour informations and modify the priority queue after an edge collapse.
         * @param vertCollapsed : the first vertex of the collapsed edge
         * @param vertCollapsedTo : the second vertex of the collapsed edge
         */
        void updateEdgeInfo(size_t vertCollapsed, size_t vertCollapsedTo);

        /**
         * @brief insert candidate edges around the given vertex into the priority queue.
         * @param vert : the id of the vertex.
         */
        void insertEdges(size_t vert);

        /**
         * @brief merge neighbour vertices are triangles of the two vertices of the collapsed edge together.
         * @param vertCollapsed : the first vertex of the collapsed edge
         * @param vertCollapsedTo : the second vertex of the collapsed edge
         */
        void mergeNeighbours(size_t vertCollapsed, size_t vertCollapsedTo);

        /**
         * @brief modify and reconstruct the triangulation after all edge collapses are finished.
         */
        void organizeOutput();

        /**
         * @brief find out the vertex that the given vertex has been collapsed to.
         * @param vert : the id of the given vertex.
         * @return the id of the vertex that the given vertex has been collapsed to.
         */
        size_t getCollapsedVert(size_t vert);

        /**
         * @brief test whether an edge collapse satisfies the link condition.
         * @param firstVert : the first vertex of the collapsed edge
         * @param secondVert : the second vertex of the collapsed edge
         * @return
         */
        bool testLinkCondition(size_t firstVert, size_t secondVert);

        /**
         * @brief build the one ring area information of an edge. The information is used to costruct the kernel region later.
         * @param firstVert : the first vertex of the collapsed edge.
         * @param secondVert : the second vertex of the collapsed edge.
         * @param boundary_faces : output boundary faces of the one ring area. This is a 3xN matrix.
         * @param related_vert_for_boundary_face : output the corresponding vertex for each boundary face of the one ring area.
         * @param innerSample : output inner shell samples that are inside the one ring area.
         * @param outerSample : output outer shell samples that are inside the one ring area.
         */
        void buildOneRingArea(size_t firstVert, size_t secondVert, matrixs_t& boundary_faces, std::vector<size_t>& related_vert_for_boundary_face,
                              std::set<size_t>& innerSample, std::set<size_t>& outerSample);

        /**
         * @brief find out the boundary faces of the one ring area of the given edge. This is used by buildOneRingArea().
         * @param firstVert : the first vertex of the collapsed edge.
         * @param secondVert : the second vertex of the collapsed edge.
         * @param boundaryFaces : output boundary faces of the one ring area. This is a 3xN matrix.
         */
        void findBoundaryFace(size_t firstVert, size_t secondVert, std::vector<std::tuple<size_t, size_t, size_t> >& boundaryFaces);

        /**
         * @brief find out the shell samples that are inside the one ring area of the given vertex. This is used by buildOneRingArea().
         * @param vert : the id of the given vertex.
         * @param innerSample : output inner shell samples that are inside the one ring area.
         * @param outerSample : output outer shell samples that are inside the one ring area.
         */
        void findShellSamples(size_t vert, std::set<size_t>& innerSample, std::set<size_t>& outerSample);

        /**
         * @brief find out the collapse position given an edge.
         *        The QEM, link condition, classification condition, normal condition are taken into consideration.
         *
         * @param vert : the first vertex of the given edge.
         * @param vertCollapseTo : the second vertex of the given edge.
         * @param position : output collapse position.
         * @param out_error : output error metric of this collapse operation.
         * @return true if found a candidate position. Otherwise, the given edge cannot be collapsed and return false.
         */
        bool findCollapsePos(size_t vert, size_t vertCollapseTo, vec3_t& position, double& out_error);

        /**
         * @brief check whether a candidate edge collapse operation satisfies the normal condition.
         * @param oneRingFaces : the boundary of the one ring area of the collapsed edge.
         * @param point : the position to collapse to.
         * @param pointType : the point type of the collapsed point.
         * @return true if satisfied, and false otherwise.
         */
        bool checkNormal(const matrixs_t& oneRingFaces, const vec3_t& point, PointType pointType);

        /**
         * @brief check whether the given point is a vertex of the given cell.
         * @param samplePoint : the position of the given point.
         * @param cellVerts : the id of the vertices of the given cell.
         * @return true if is, and false otherwise.
         */
        bool isVertOfCell(const vec3_t& samplePoint, const std::vector<size_t>& cellVerts);

        /**
         * @brief check whether the given cell has already been collapsed.
         * @param cellVerts : the id of the vertices of the given cell.
         * @return true if is, and false otherwise.
         */
        bool isCollapsedCell(const std::vector<size_t>& cellVerts);

    private:
        TriangulatedShell& mTriangulation;
        const Shell& mShell;
        Type mType;
        bool mIsHalfEdge = true;
        double mSampleRadius;

        std::vector<matrixr_t> mQ;  /**< store the Q matrices for each vertex*/
        std::vector<std::set<size_t> > mNeighbours;   /**< store the neighbour vertices for each vertex */
        std::vector<std::set<size_t> > mCollapseableNeighbours;   /**< store the collapseable neighbour vertices for each vertex */
        std::vector<std::set<std::pair<size_t, size_t> > > mSameShellNeighbourFaces;  /**< store the neighbour faces which are on same shell with the vertex. */
        std::vector<std::set<size_t> > mNeighbourCells; /**< store neighbour cells for each vertex */

        //for half edge collapse
        std::vector<size_t> mCollapseTo;    //record what vertex the vertices has collapsed to.
        std::vector<std::vector<std::shared_ptr<EdgeInfo> > > mRelatedEdgeInfo;    //recording the edge infos related to the vertices

        std::mutex mtx; // used for synchronizing insertEdges() between threads.

    };
}

#endif
