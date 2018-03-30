#ifndef WKY_SBV_FIELD_COMPUTER_H
#define WKY_SBV_FIELD_COMPUTER_H

#include "Common.h"

namespace SBV
{
    // information of triangle for integrating Green Function
    struct Triangle
    {
        vec3_t position; // barycenter of the triangle
        vec3_t normal;   // normal of the triangle
        double value = 0; // boundary value on the triangle
        double derivative = 0; // boundary derivative on the triangle
        double size = 0;   // size of the triangle.
        mat3x3_t tri;      // the triangle
        mat4x4_t transform;  //matrix for transforming to local coordinate described in Graglia's paper
        mat4x4_t invTransform; // inverse matrix of the above transform matrix
    };

    class FieldComputerImpl;

    class FieldComputer{
    public:
        FieldComputer() = default;
        ~FieldComputer();

        void init(const std::vector<Triangle>& samples);
        double getFieldValue(const vec3_t& x);

    private:
        FieldComputerImpl* impl = nullptr;
    };
}

#endif
