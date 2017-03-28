#ifndef WKY_SIMPLIFIER_H
#define WKY_SIMPLIFIER_H

#include <string>
#include "Common.h"
#include "Refinement.h"

namespace SBV {
    class Simplifier{
    public:
        Simplifier(Curve& mesh);

        void simplify();

        void setOutputDirectory(const std::string& outputDir)
        {
            this->mOutputDirectory = outputDir;
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
            matrixr_t vertices;
            matrixs_t lines;
        };

    private:
        void genDefaultParams();
        void generateShells();
        void sample(const matrixr_t& vertices, const matrixs_t& triangles, std::vector<matrixr_t>& output_samples);
        void refine();
        void collapseBoundary();
        void mutualTessellate();

    private:
        std::string mOutputDirectory;
        Curve& mSourceMesh;
        double mMaxDistance = -1;
        double mSampleRadius = -1;
        double mAlpha = 0.2;
        bool mNeedGenTempResult = false;

        matrixr_t mInnerShell;
        matrixr_t mOuterShell;

        TriangulatedShell mTriangulation;
        ZeroSet mZeroSet;
    };
}

#endif
