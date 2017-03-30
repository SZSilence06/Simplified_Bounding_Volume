#ifndef WKY_QUADTREE_H
#define WKY_QUADTREE_H

#include <Geometry/Common.h>
#include <Common.h>

namespace  WKYLIB
{
    namespace Geometry
    {
        class QuadTree
        {
        public:
            QuadTree(double xMin, double xMax, double yMin, double yMax);
            QuadTree(double xMin, double xMax, double yMin, double yMax, int maxDepth);
            ~QuadTree();

            void insert(const matrixr_t& point);
            void getNearestPoint(const matrixr_t& point, matrixr_t& nearest);

        private:
            struct QuadTreeNode
            {
                enum Type
                {
                    NODE_EMPTY,
                    NODE_LEAF,
                    NODE_POINTER
                };

                matrixr_t point;
                QuadTreeNode* lt = nullptr;
                QuadTreeNode* lb = nullptr;
                QuadTreeNode* rt = nullptr;
                QuadTreeNode* rb = nullptr;
                QuadTreeNode* parent = nullptr;
                double xMax;
                double xMin;
                double yMax;
                double yMin;
                int depth = 0;
                Type type = NODE_EMPTY;

                QuadTreeNode();
                QuadTreeNode(double xMin, double xMax, double yMin, double yMax, QuadTreeNode* parent);
            };

        private:
            void insertToNode(const matrixr_t& point, QuadTreeNode* node);
            void split(QuadTreeNode* node);
            QuadTreeNode* getNodeForInsert(const matrixr_t& point, QuadTreeNode* parent);

        private:
            QuadTreeNode* mRoot = nullptr;

            double mXMax;
            double mXMin;
            double mYMax;
            double mYMin;
            double mMaxDepth;
        };
    }
}

#endif
