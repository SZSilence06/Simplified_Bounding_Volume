#ifndef WKY_SBV_FIELD_COMPUTER_H
#define WKY_SBV_FIELD_COMPUTER_H

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
        mat4x4_t tri;      //indicating the triangle which the sample point lies in.
        mat4x4_t transform;  //matrix for transforming to local
        mat4x4_t invTransform;
    };

    class FieldComputerImpl;

    class FieldComputer{
    public:
        FieldComputer() = default;
        ~FieldComputer();

        void init(const std::vector<SamplePoint>& samples);
        double getFieldValue(const vec3_t& x);

    private:
        FieldComputerImpl* impl = nullptr;
    };
}

#endif
