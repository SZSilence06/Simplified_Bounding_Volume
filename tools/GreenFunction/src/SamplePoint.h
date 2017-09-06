#ifndef GF_SAMPLE_POINT_H
#define GF_SAMPLE_POINT_H

#include "Common.h"

struct SamplePoint
{
    vec2_t position;
    matrixr_t normal;
    double color = 0;
    double derivative = 0;
};

#endif
