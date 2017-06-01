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

#ifdef VER_2D
        Simplifier(Curve& mesh);
#else
        Simplifier(Mesh& mesh);
#endif

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
            std::vector<Point> vertices;
            std::vector<Eigen::Vector2i> lines;
        };

    private:
        void genDefaultParams();
        void generateShells();
#ifdef VER_2D
        void sample(const std::vector<Point> &vertices, const std::vector<Eigen::Vector2i> &lines,
                    std::vector<Point> &output_samples);
#else
        void sample(const std::vector<Point> &vertices, const std::vector<Eigen::Vector3i> &triangles,
                    std::vector<Point> &output_samples);
#endif
        void refine();
        void collapseBoundary();
        void mutualTessellate();
        void collapseZeroSet();
        void writeSummary();
        void logTimer(const std::string& prefix, const DebugTimer& timer);

    private:
        std::string mOutputDirectory;
#ifdef VER_2D
        Curve& mSourceMesh;
#else
        Mesh& mSourceMesh;
#endif

        double mMaxDistance = -1;
        double mSampleRadius = -1;
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
