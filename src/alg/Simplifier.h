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

        void setGenTempResult(bool value)
        {
            this->mNeedGenTempResult = value;
        }

    public:
        using TriangulatedShell = Refinement::TriangulatedShell;

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
        void buildZeroSet();
        size_t getZeroPointIndex(size_t firstVertex, size_t secondVertex, std::vector<std::pair<size_t, size_t> >& existingVertPairs);

    private:
        std::string mOutputDirectory;
        Curve& mSourceMesh;
        double mMaxDistance = -1;
        double mSampleRadius = -1;
        bool mNeedGenTempResult = false;

        matrixr_t mInnerShell;
        matrixr_t mOuterShell;

        TriangulatedShell mTriangulation;
        ZeroSet mZeroSet;
    };
}

#endif
