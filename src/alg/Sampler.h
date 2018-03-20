#ifndef WKY_SAMPLER_H
#define WKY_SAMPLER_H

#include "Common.h"

namespace SBV
{
    /**
     * @brief Static class to perform sampling algorithms.
     */
    class Sampler
    {
    public:
        /**
         * @brief The poisson disk sampling.
         * @param vertices : vertex matrix of the sampled mesh.
         * @param triangles : triangle matrix of the sampled mesh.
         * @param radius : sample radius.
         * @param samples : output samples.
         * @param sample_normals : output sample normals.
         */
        static void poissonDisk(const matrixr_t& vertices, const matrixs_t& triangles, double radius,
                            matrixr_t& samples, matrixr_t& sample_normals);
    };
}

#endif
