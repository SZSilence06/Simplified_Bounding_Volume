#include "QuadTree.h"

namespace WKYLIB
{
    namespace Geometry
    {
        QuadTree::QuadTree(double xMin, double xMax, double yMin, double yMax)
            : QuadTree(xMin, xMax, yMin, yMax, 10000)
        {

        }

        QuadTree::QuadTree(double xMin, double xMax, double yMin, double yMax, int maxDepth)
            : mXMin(xMin),
              mXMax(xMax),
              mYMin(yMin),
              mYMax(yMax),
              mMaxDepth(maxDepth)
        {

        }

        QuadTree::~QuadTree()
        {
            if(mRoot)
            {
                delete mRoot;
            }
        }

        void QuadTree::insert(const matrixr_t &point)
        {
            if(point[0] < mXMin || point[0] > mXMax || point[1] < mYMin || point[1] < mYMax)
            {
                throw std::invalid_argument("the point to insert is beyond the bounding box of quadtree.");
            }
            insertToNode(point, mRoot);
        }

        void QuadTree::getNearestPoint(const matrixr_t &point, matrixr_t &nearest)
        {

        }

        void QuadTree::insertToNode(const matrixr_t &point, QuadTreeNode *node)
        {
            switch(node->type)
            {
            case QuadTreeNode::NODE_EMPTY:
                node->point = point;
                node->type = QuadTreeNode::NODE_LEAF;
                break;
            case QuadTreeNode::NODE_LEAF:
                split(node);
                insertToNode(point, node);
                break;
            case QuadTreeNode::NODE_POINTER:
                insertToNode(point, getNodeForInsert(point, node));
                break;
            }
        }

        void QuadTree::split(QuadTreeNode *node)
        {
            assert(node->type == QuadTreeNode::NODE_LEAF);

            double xMid = (node->xMin + node->xMax) / 2;
            double yMid = (node->yMin + node->yMax) / 2;

            node->lb = new QuadTreeNode(node->xMin, xMid, node->yMin, yMid, node);
            node->lt = new QuadTreeNode(node->xMin, xMid, yMid, node->yMax, node);
            node->rt = new QuadTreeNode(xMid, node->xMax, node->yMin, yMid, node);
            node->rb = new QuadTreeNode(xMid, node->xMax, yMid, node->yMax, node);
            node->type = QuadTreeNode::NODE_POINTER;

            insertToNode(node->point, node);
        }

        QuadTree::QuadTreeNode* QuadTree::getNodeForInsert(const matrixr_t &point, QuadTreeNode *parent)
        {
            assert(parent->type == QuadTreeNode::NODE_POINTER);

            double xMid = (parent->xMin + parent->xMax) / 2;
            double yMid = (parent->yMin + parent->yMax) / 2;
            bool isL = point[0] < xMid;
            bool isB = point[1] < yMid;
            if(isL)
            {
                return isB ? parent->lb : parent->lt;
            }
            return isB ? parent->rb : parent->rt;
        }

        /////////////////////////////////////////////
        QuadTree::QuadTreeNode::QuadTreeNode()
        {

        }

        QuadTree::QuadTreeNode::QuadTreeNode(double xMin, double xMax, double yMin, double yMax, QuadTreeNode *parent)
            : xMin(xMin),
              xMax(xMax),
              yMin(yMin),
              yMax(yMax),
              parent(parent),
              depth(parent->depth + 1)
        {

        }
    }
}
