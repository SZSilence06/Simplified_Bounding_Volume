#ifndef GF_CURVE_H
#define GF_CURVE_H

#include "Common.h"

struct Curve
{
    matrixr_t points;
    matrixs_t lines;

    void adjust();

private:
    void buildAdjacency(std::vector<std::vector<int>>& adjacency);
};



#endif
