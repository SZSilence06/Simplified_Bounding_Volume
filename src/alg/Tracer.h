#ifndef WKY_SBV_TRACER_H
#define WKY_SBV_TRACER_H

#include "Common.h"
#include "FieldComputer.h"

namespace SBV
{
    class Tracer{
    public:
        Tracer(const std::vector<SamplePoint>& samples) : mSamples(samples) {}

        void tracePoints(const matrixr_t& points, const matrixr_t& normals, double length, matrixr_t& traceResult);

    private:
        const std::vector<SamplePoint>& mSamples;
    };
}

#endif
