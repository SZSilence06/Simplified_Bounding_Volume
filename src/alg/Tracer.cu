#include "Tracer.h"
#include <zjucad/matrix/matrix.h>

using namespace zjucad::matrix;

namespace SBV {
    template<class T>
    __device__ static inline T my_min(const T& a, const T& b) {
        return a < b ? a : b;
    }

    __device__ static void integrateOverTriangle(const vec3_t& x, const SamplePoint &point, double& I1, vec3_t& Igrad)
    {
        mat4x4_t triangle;
        triangle(colon(0, 2), colon(0, 2)) = point.tri;
        triangle(3, colon(0, 2)) = ones<double>(1, 3);
        mat4x4_t localTriangle = point.transform * triangle;

        double l3 = localTriangle(0, 1);
        double u3 = localTriangle(0, 2);
        double v3 = localTriangle(1, 2);

        vec4_t tempX;
        tempX(colon(0, 2), colon()) = x;
        tempX[3] = 1;
        vec4_t localX = point.transform * tempX;
        double u0 = localX[0];
        double v0 = localX[1];
        double w0 = localX[2];

        // edge lengths
        double l1 = sqrt((l3-u3) * (l3-u3) + v3*v3);
        double l2 = sqrt(u3*u3 + v3*v3);

        // threshold for small numbers
        double threshold = 1e-6 * my_min(my_min(l1,l2), l3);
        if(fabs(w0) < threshold)
            w0 = 0;

        // eq (3)
        vec3_t sminus, splus;
        sminus[0] = -((l3-u3)*(l3-u0)+v3*v0)/l1;
        sminus[1] = -(u3*(u3-u0)+v3*(v3-v0))/l2;
        sminus[2] = -u0;
        splus[0] = ((u3-l3)*(u3-u0)+v3*(v3-v0))/l1;
        splus[1] = (u3*u0+v3*v0)/l2;
        splus[2] = l3-u0;

        // eq (4)
        vec3_t t0;
        t0[0] = ((u3-l3)*v0+v3*(l3-u0))/l1;
        t0[1] = (v3*u0-u3*v0)/l2;
        t0[2] = v0;

        // eq (5)
        vec3_t tplus, tminus;
        tplus[0] = sqrt((u3-u0)*(u3-u0) + (v3-v0)*(v3-v0));
        tplus[1] = sqrt(u0*u0 + v0*v0);
        tplus[2] = sqrt((l3-u0)*(l3-u0) + v0*v0);
        tminus[0] = tplus[2];
        tminus[1] = tplus[0];
        tminus[2] = tplus[1];

        // line 1, pp. 1450
        vec3_t R0;
        for(int i = 0; i < 3; i++)
            R0[i] = sqrt(t0[i]*t0[i] + w0*w0);

        //line 2, pp. 1450
        vec3_t Rplus, Rminus;
        for(int i = 0; i < 3; i++)
        {
            Rplus[i] = sqrt(tplus[i]*tplus[i] + w0*w0);
            Rminus[i] = sqrt(tminus[i]*tminus[i] + w0*w0);
        }

        // eq (11)
        vec3_t f2;
        for(int i = 0; i < 3; i++)
        {
            double temp;
            if(w0 == 0)
            {
                if(fabs(t0[i]) < threshold)
                    temp = fabs(log(splus[i]) / sminus[i]);
                else
                    temp = (tplus[i]+splus[i]) / (tminus[i]+sminus[i]);
            }
            else
            {
                 temp = (Rplus[i]+splus[i]) / (Rminus[i]+sminus[i]);
            }
            f2[i] = log(temp);
            //fix value for points on the triangle corners
            if(f2[i] != f2[i])  //nan
                f2[i] = 0;
        }


        // eq (13) and eq (14)
        vec3_t beta;
        double betaSum;
        if(w0 == 0)
        {
            for(int i = 0; i < 3; i++)
            {
                if(fabs(t0[i]) < threshold)
                    beta[i] = 0;
                else
                    beta[i] = atan(splus[i] / t0[i]) - atan(sminus[i] / t0[i]);
            }
        }
        else
        {
            for(int i = 0; i < 3; i++)
                beta[i] = atan((t0[i]*splus[i]) / (R0[i]*R0[i] + Rplus[i]*fabs(w0))) - atan((t0[i]*sminus[i]) / (R0[i]*R0[i] + Rminus[i]*fabs(w0)));
        }
        betaSum = beta[0] + beta[1] + beta[2];


        // eq (19), integral of kernel 1/R
        I1 = 0;
        for(int i = 0; i < 3; i++)
            I1 += t0[i]*f2[i];
        I1 -= fabs(w0) * betaSum;

        // normals of the triangle edges, Fig. 1(b)
        vec3_t m[3];
        m[0][0] = v3;
        m[0][1] = l3 - u3;
        m[0][2] = 0;
        m[1][0] = -v3;
        m[1][1] = u3;
        m[1][2] = 0;
        m[2][0] = 0;
        m[2][1] = -l3;
        m[2][2] = 0;
        for(int i = 0; i < 3; i++)
            m[i] /= norm(m[i]);

        // eq (34), integral of kernel grad(1/R)
        Igrad = zeros<double>(3, 1);
        for(int i = 0; i < 3; i++)
            Igrad -= m[i] * f2[i];
        vec3_t w = zeros<double>(3, 1);
        w[2] = 1;
        if(w0 >= 0)
            Igrad -= betaSum * w;
        else if(w0 < 0)
            Igrad += betaSum * w;

        //transform back to world space
        vec4_t IgradTemp;
        IgradTemp(colon(0, 2), 0) = Igrad;
        IgradTemp[3] = 0;
        vec4_t IgradGlob = point.invTransform * IgradTemp;
        Igrad = IgradGlob(colon(0, 2), 0);
    }

    __device__ static double kernel(const vec3_t& x, const SamplePoint& sample)
    {
        double I1;
        vec3_t Igrad;
        integrateOverTriangle(x, sample, I1, Igrad);
        //double result = (sample.value * dot(Igrad, sample.normal) - sample.derivative * I1) / (-4 * PI);
        double result = I1 * sample.derivative;
        return result;
    }

    __device__ static double getFieldValue(const SamplePoint* samples, int* sampleCount, const vec3_t &x)
    {
        double result = 0;
        for(int i = 0; i < *sampleCount; i++)
        {
            result += kernel(x, samples[i]);
        }
        return result;
    }

    __device__ vec3_t getGradient(const SamplePoint* samples, int* sampleCount, const vec3_t &x)
    {
        const double STEP = 0.0001;
        double a = getFieldValue(samples, sampleCount, x);

        vec3_t xx = x;
        xx[0] += STEP;
        double b = getFieldValue(samples, sampleCount, xx);

        vec3_t xy = x;
        xy[1] += STEP;
        double c = getFieldValue(samples, sampleCount, xy);

        vec3_t xz = x;
        xz[2] += STEP;
        double d = getFieldValue(samples, sampleCount, xz);

        vec3_t result;
        result[0] = (b - a) / STEP;
        result[1] = (c - a) / STEP;
        result[2] = (d - a) / STEP;
        return result;
    }

    __global__ void kernel_trace(const SamplePoint* samples, int* sampleCount, vec3_t* points, vec3_t* normals, int* pointCount, double* distance)
    {
        const double STEP = 0.01;
        double t = 0;
        const int INDEX = blockIdx.x;
        vec3_t result = points[INDEX];
        vec3_t n = normals[INDEX];
        while(t < *distance)
        {
            if(t == 0)
            {
                result += STEP * n;
            }
            else
            {
                vec3_t grad = getGradient(samples, sampleCount, result);
                double norm_grad = norm(grad);
                grad /= norm_grad;
                result -= STEP * grad;
            }
            t += STEP;
        }
        points[INDEX] = result;
    }
}

