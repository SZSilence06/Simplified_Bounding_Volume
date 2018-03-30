#ifndef WKY_SIMPLIFIER_H
#define WKY_SIMPLIFIER_H

#include <string>
#include "Common.h"
#include "Refinement.h"
#include "Shell.h"
#include "Logger.h"
#include <wkylib/debug_util.h>

namespace SBV {
     /**
     * @brief The Simplifier class
     *
     * This class is used to perform the reconstruction and simplification of the mesh from input sample points.
     * The sample points are organized as an inner shell and an outer shell.
     */
    class Simplifier{
    public:
        using DebugTimer = WKYLIB::DebugTimer;

        /**
         * @brief Simplifier constructor
         * @param shell : the input sample points.
         */
        Simplifier(Shell& shell);

        /**
         * @brief generate the result mesh.
         */
        void simplify();

        /**
         * @brief set the output directory of the result.
         * @param outputDir : output directory.
         */
        void setOutputDirectory(const std::string& outputDir)
        {
            this->mOutputDirectory = outputDir;
            Logger::getInstance().setFile(outputDir + "/log.txt");
        }

        /**
         * @brief set the alpha parameter. For more information about this parameter, see the paper.
         * @param alpha : the alpha parameter.
         */
        void setAlpha(double alpha)
        {
            this->mAlpha = alpha;
        }

        /**
         * @brief set whether to output intermediate results.
         * @param value : true to output, false otherwise.
         */
        void setGenTempResult(bool value)
        {
            this->mNeedGenTempResult = value;
        }   

    public:
        struct ZeroSet
        {
            matrixr_t vertices;
            matrixs_t lines;
        };

    private:
        void genDefaultParams();
        void refine();
        void collapseBoundary();
        void mutualTessellate();
        void collapseZeroSet();
        void writeSummary();
        void logTimer(const std::string& prefix, const DebugTimer& timer);

    private:
        std::string mOutputDirectory;
        double mAlpha = 0.2;
        bool mNeedGenTempResult = false;

        DebugTimer mTimerRefine = DebugTimer("refinement");
        DebugTimer mTimerBoundaryHalfEdge = DebugTimer("Boundary Collapse(Half Edge)");
        DebugTimer mTimerBoundaryGeneral = DebugTimer("Boundary Collapse(General)");
        DebugTimer mTimerZeroSetHalfEdge = DebugTimer("Zero Set Collapse(Half Edge)");
        DebugTimer mTimerZeroSetGeneral = DebugTimer("Zero Set Collapse(General)");
        DebugTimer mTimerMutualTessellation = DebugTimer("Mutual Tesselation");
        DebugTimer mTimerSimplify = DebugTimer("Simplification");

        Shell mShell;

        TriangulatedShell mTriangulation;
        ZeroSet mZeroSet;
    };
}

#endif
