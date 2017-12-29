#ifndef WKY_SIMPLIFIER_H
#define WKY_SIMPLIFIER_H

#include <string>
#include "Common.h"
#include "Refinement.h"
#include "Shell.h"
#include "Logger.h"
#include <wkylib/debug_util.h>

namespace SBV {
    class Simplifier{
    public:
        using DebugTimer = WKYLIB::DebugTimer;

        Simplifier(Shell& shell);

        void simplify();

        void setOutputDirectory(const std::string& outputDir)
        {
            this->mOutputDirectory = outputDir;
            Logger::getInstance().setFile(outputDir + "/log.txt");
        }

        void setMaxDistance(double maxDist)
        {
            this->mMaxDistance = maxDist;
        }

        void setSampleRadius(double sampleRadius)
        {
            this->mSampleRadius = sampleRadius;
        }

        void setSampleCount(int sampleCount)
        {
            this->mSampleCount = sampleCount;
        }

        void setAlpha(double alpha)
        {
            this->mAlpha = alpha;
        }

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
        void sample(const matrixr_t& vertices, const matrixs_t& triangles, matrixr_t& output_samples, matrixr_t& normals);
        void refine();
        void collapseBoundary();
        void mutualTessellate();
        void collapseZeroSet();
        void writeSummary();
        void logTimer(const std::string& prefix, const DebugTimer& timer);

    private:
        std::string mOutputDirectory;
        double mMaxDistance = -1;
        double mSampleRadius = -1;
        double mAlpha = 0.2;
        bool mNeedGenTempResult = false;
        int mSampleCount = -1;

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
