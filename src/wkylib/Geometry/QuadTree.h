#ifndef WKY_GEOMETRY_QUADTREE_H
#define WKY_GEOMETRY_QUADTREE_H

#include "Common.h"
#include "../Common.h"

namespace  WKYLIB
{
    namespace Geometry
    {
        class QuadTreeIterator;

        class QuadTree
        {
        public:
            QuadTree(double xMin, double xMax, double yMin, double yMax);
            QuadTree(double xMin, double xMax, double yMin, double yMax, int maxDepth);
            ~QuadTree();

            void insert(const matrixr_t& point);
            void getNearestPoint(const matrixr_t& point, matrixr_t& nearest);

            using Iterator = QuadTreeIterator;

        protected:
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
                ~QuadTreeNode();
            };

        protected:
            void insertToNode(const matrixr_t& point, QuadTreeNode* node);
            void split(QuadTreeNode* node);
            QuadTreeNode* getNodeForInsert(const matrixr_t& point, QuadTreeNode* parent);

        protected:
            QuadTreeNode* mRoot = nullptr;

            double mXMax;
            double mXMin;
            double mYMax;
            double mYMin;
            double mMaxDepth;

            friend class QuadTreeIterator;
        };

        class QuadTreeIterator
        {
        public:
            QuadTreeIterator(const QuadTree& tree);
            ~QuadTreeIterator() = default;

            bool next(matrixr_t& point);

        private:
            enum ParentState
            {
                NODE_LB,
                NODE_LT,
                NODE_RB,
                NODE_RT,
                NODE_NULL
            };

            using QuadTreeNode = QuadTree::QuadTreeNode;

        private:
            bool findNext(QuadTreeNode* node, matrixr_t& point);
            QuadTreeNode* getNextNode(QuadTreeNode* node);
            QuadTreeNode* getNextLeafNode(QuadTreeNode* node);
            ParentState getParentState(QuadTreeNode* node);

        private:
            const QuadTree& mTree;
            QuadTreeNode* mNode;

            friend class QuadTree;
        };
    }
}

#endif
