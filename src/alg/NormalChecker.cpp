#include "NormalChecker.h"
#include "Shell.h"
#include "BaryComputer.h"

using namespace zjucad::matrix;

namespace SBV
{
    //private tool functions
    static bool checkClassification(const BaryComputer& baryComputer, PointType type_v1, PointType type_v2, PointType type_v3, PointType type_v4,
                                    const vec3_t& point, bool isOuter)
    {
        vec4_t bary;
        baryComputer(point, bary);
        double fvalue = TriangulatedShell::getFValue(type_v1) * bary[0]
                + TriangulatedShell::getFValue(type_v2) * bary[1]
                + TriangulatedShell::getFValue(type_v3) * bary[2]
                + TriangulatedShell::getFValue(type_v4) * bary[3];
        bool result;
        if(isOuter)
        {
            result = fvalue > 0;
        }
        else
        {
            result = fvalue < 0;
        }
        return result;
    }

    /////////////////////////////////////////////////////////////////////////////////
    bool NormalChecker::check(const matrixr_t &cell, PointType type_v1, PointType type_v2, PointType type_v3, PointType type_v4,
                              const Shell &shell)
    {
        BaryComputer baryComputer(cell);

        matrixr_t center = (cell(colon(), 0) + cell(colon(), 1) + cell(colon(), 2) + cell(colon(), 3)) / 4.0;
        constexpr double k = cbrt(0.7);

        matrixr_t newVerts[4];
        for(int i = 0; i < 4; i++)
        {
            newVerts[i] = k * cell(colon(), i) + (1 - k) * center;
        }

        //check inner verts classification
        const KdTreeWrap& innerTree = shell.getInnerTree();
        size_t nearVerts[4];
        for(int i = 0; i < 4; i++)
        {
            nearVerts[i] = innerTree.getNearestPoint(newVerts[i]);
        }

        for(int i = 0; i < 4; i++)
        {
            if(checkClassification(baryComputer, type_v1, type_v2, type_v3, type_v4,
                                   shell.mInnerShell(colon(), nearVerts[i]), false) == false)
            {
                return false;
            }
        }

        //check outer verts classification
        const KdTreeWrap& outerTree = shell.getOuterTree();
        for(int i = 0; i < 4; i++)
        {
            nearVerts[i] = outerTree.getNearestPoint(newVerts[i]);
        }

        for(int i = 0; i < 4; i++)
        {
            if(checkClassification(baryComputer, type_v1, type_v2, type_v3, type_v4,
                                   shell.mOuterShell(colon(), nearVerts[i]), true) == false)
            {
                return false;
            }
        }

        return true;
    }
}
