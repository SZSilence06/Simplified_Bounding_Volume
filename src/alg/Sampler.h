#ifndef WKY_SAMPLER_H
#define WKY_SAMPLER_H

#include "Common.h"

namespace SBV
{
    class Sampler
    {
    public:
        static void poisson(const matrixr_t& vertices, const matrixs_t& triangles, double radius,
                            matrixr_t& samples, matrixr_t& sample_normals);
    };
}

#endif
