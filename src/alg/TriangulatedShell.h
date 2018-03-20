#ifndef WKY_TRIANGULATED_SHELL_H
#define WKY_TRIANGULATED_SHELL_H

#include "Common.h"
#include <set>

namespace SBV
{
    enum PointType
    {
        POINT_BOUNDING_BOX,
        POINT_INNER,
        POINT_OUTER,
        POINT_ZERO,
        POINT_UNKNOWN
    };

    /**
     * @brief The TriangulatedShell class
     *
     * This class describes the delauny triangulation.
     */
    class TriangulatedShell
    {
    public:
        struct ZeroSet
        {
            matrixr_t vertices;
            matrixs_t triangles;
            std::vector<std::pair<size_t, size_t> >  vertPairs;   //recording the corresponding vert pairs of the zero set vertices
            std::vector<size_t> faceTetras;                        //recording the corresponding tetras of the zero set triangles
        };

    public:
        matrixr_t vertices;
        matrixs_t cells;
        std::vector<PointType> vertType;

    public:
        /**
         * @brief get the F function value of the given point type.
         * @param pointType : the given point type.
         * @return the F function value.
         */
        static double getFValue(PointType pointType);

        /**
         * @brief get the F function value of the given vertex.
         * @param vert : the id of the given vertex.
         * @return the F function value.
         */
        double getFValue(size_t vert) const;

        /**
         * @brief get the sign of the F function value of the given vertex.
         * @param vert : the id of the given vertex.
         * @return 1 if >0, 0 if =0 and -1 if < 0.
         */
        double getSign(size_t vert) const;

        /**
         * @brief returns the zero set.
         */
        inline const ZeroSet& getZeroSet() const { return mZeroSet; }

        /**
         * @brief extract zero set from the triangulation.
         */
        void buildZeroSet();

        /**
         * @brief do the mutual tessellation.
         */
        void mutualTessellate();

        /**
         * @brief write the inner shell and outer shell triangles of the triangulation.
         * @param path : the path to the output directory.
         */
        void outputBoundaryFaces(const std::string& path);

    private:
        /**
         * @brief This struct describes a face on the zero set. This is used by buildZeroSetExisting().
         */
        struct ZeroFace{
            std::set<size_t> verts; // records the id of the vertices of this zero set face.
            size_t tetra;          // record the corresponding tetrahedron of this zero set face.

            bool operator <(const ZeroFace& face) const { return this->verts < face.verts; } // compare function used by std::set<ZeroFace>
        };

        /**
         * @brief get the actual zero set point id of the given edge. If the zero point is not contained in the current zero set,
         *        insert it.
         *
         *        This is because an edge is shared by multiple tetrahedrons, so same zero point may be queried for multiple times.
         *        If each time when we extract zero point from an edge we always create a new zero set point and insert to the zero set,
         *        then the zero set will contain many duplicated vertices.
         *        So we maintain a set and uses this function to query whether the zero point is already inserted.
         *        If so, return the id. Otherwise, insert the point and then return the id.
         *
         *
         * @param firstVertex : the id of the first vertex of the edge to extract zero set point.
         * @param secondVertex : the id of the second vertex of the edge to extract zero set point.
         * @return The actual zero set point id of the given edge.
         */
        size_t getZeroPointIndex(size_t firstVertex, size_t secondVertex);

        /**
         * @brief buildZeroSetExisting
         *
         * Extract the zero set if the triangulation already contains zero set points. This happens after the mutual tesselation.
         */
        void buildZeroSetExisting();

        /**
         * @brief try to add a zero set face into the zeroFaces if not existed.
         * @param currentTetra : the id of the tetrahedron corresponding to the three given vertices.
         * @param zeroVert1 : the id of the first given vertex
         * @param zeroVert2 : the id of the second given vertex
         * @param zeroVert3 : the id of the third given vertex
         * @param zeroFaces : the set of zero set faces to insert
         */
        void tryAddZeroFace(size_t currentTetra, size_t zeroVert1, size_t zeroVert2, size_t zeroVert3, std::set<ZeroFace>& zeroFaces);

    private:
        ZeroSet mZeroSet;
        bool hasZeroSet = false;
    };
}

#endif
