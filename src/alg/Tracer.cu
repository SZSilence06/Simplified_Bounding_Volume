#include "Tracer.h"
#include <zjucad/matrix/matrix.h>
#include <wkylib/Cuda/CudaPointer.h>
#include <wkylib/Cuda/CudaVector.h>
#include <eigen3/Eigen/Dense>
#include <iostream>

using namespace zjucad::matrix;
using namespace WKYLIB::Cuda;

namespace SBV {
    static vec3_t getGradient(FieldComputer* field, const vec3_t &x)
    {
        const double STEP = 0.0001;
        double a = field->getFieldValue(x);

        vec3_t xx = x;
        xx[0] += STEP;
        double b = field->getFieldValue(xx);

        vec3_t xy = x;
        xy[1] += STEP;
        double c = field->getFieldValue(xy);

        vec3_t xz = x;
        xz[2] += STEP;
        double d = field->getFieldValue(xz);

        vec3_t result;
        result[0] = (b - a) / STEP;
        result[1] = (c - a) / STEP;
        result[2] = (d - a) / STEP;
        return result;
    }

    void Tracer::tracePoints(const matrixr_t &points, const matrixr_t &normals, double length, matrixr_t &traceResult)
    {
        const double STEP = length / 10;

        traceResult.resize(3, points.size(2));

        for(int i = 0; i < points.size(2); i++) {
            std::cout << "Tracing point " << i << " / " << points.size(2) << std::endl;
            vec3_t result = points(colon(), i);
            vec3_t n = normals(colon(), i);
            double t = STEP;
            result += STEP * n;
            while(t < length)
            {
                vec3_t grad = getGradient(mField.get(), result);
                double norm_grad = norm(grad);
                grad /= norm_grad;
                result -= STEP * grad;
                t += STEP;
            }
            traceResult(colon(), i) = result;
        }
    }
}

