#ifndef WKY_SBV_TRACER_H
#define WKY_SBV_TRACER_H

#include "Common.h"

namespace SBV
{
    //for computing Green Function
    struct SamplePoint
    {
        vec3_t position;
        vec3_t normal;
        double value = 0;
        double derivative = 0;
        double size = 0;   //indicating the size of the triangle which the sample point lies in.
        mat3x3_t tri;      //indicating the triangle which the sample point lies in.
        matrixr_t transform;  //matrix for transforming to local
        matrixr_t invTransform;
    };

    class Tracer{
    public:
        Tracer(const std::vector<SamplePoint>& samples) : mSamples(samples) {}

        void tracePoints(const matrixr_t& points, const matrixr_t& normals, double length, matrixr_t& traceResult);

    private:
        const std::vector<SamplePoint>& mSamples;
    };
}

#endif
