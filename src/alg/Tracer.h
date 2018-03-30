#ifndef WKY_SBV_TRACER_H
#define WKY_SBV_TRACER_H

#include "Common.h"
#include "FieldComputer.h"

namespace SBV
{
     /**
     * @brief The Tracer class
     *
     *        This class is used for tracing points along the gradient line.
     *
     *        There are two possible parallel strategies:
     *
     *        1. Uses cuda to compute the field value integration parallelly, and trace points sequencially. This is the current strategy.
     *        2. Uses cuda to trace points parallelly, and compute the field with fmm acceleration.
     */
    class Tracer{
    public:
        Tracer(const std::vector<Triangle>& samples) : mSamples(samples) {
            mField = std::unique_ptr<FieldComputer>(new FieldComputer());
            mField->init(samples);
        }

        void tracePoints(const matrixr_t& points, const matrixr_t& normals, double length, matrixr_t& traceResult);

    private:
        const std::vector<Triangle>& mSamples;

        std::unique_ptr<FieldComputer> mField;
    };
}

#endif
