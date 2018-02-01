#include "Util.h"

namespace SBV{
    template<class T>
    __host__ __device__ static inline T my_min(const T& a, const T& b) {
        return a < b ? a : b;
    }

    __host__ __device__ static void viewTransform(const Eigen::Vector3d &eye, Eigen::Vector3d ux, Eigen::Vector3d uz, Eigen::Matrix4d &output)
    {
        uz /= uz.norm();
        ux /= ux.norm();
        const Eigen::Vector3d uy = uz.cross(ux);

        output(0, 0) = ux[0];
        output(0, 1) = ux[1];
        output(0, 2) = ux[2];
        output(0, 3) = -eye.dot(ux);
        output(1, 0) = uy[0];
        output(1, 1) = uy[1];
        output(1, 2) = uy[2];
        output(1, 3) = -eye.dot(uy);
        output(2, 0) = uz[0];
        output(2, 1) = uz[1];
        output(2, 2) = uz[2];
        output(2, 3) = -eye.dot(uz);
        output(3, 0) = 0;
        output(3, 1) = 0;
        output(3, 2) = 0;
        output(3, 3) = 1;
    }

    __host__ __device__ static void localTransform(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, Eigen::Matrix4d &output)
    {
        Eigen::Vector3d ux = b - a;
        Eigen::Vector3d uz = ux.cross(c - a);
        viewTransform(a, ux, uz, output);
    }

    //closed-form calculation according to Graglia 1993.
    __host__ __device__ double Util::GPU_integrateOverTriangle(const Eigen::Vector3d& x, const Eigen::Matrix3d &triangle)
    {
        Eigen::Matrix4d triangle_4;
        triangle_4.block<3, 3>(0, 0) = triangle;
        triangle_4(3, 0) = triangle_4(3, 1) = triangle_4(3, 2) = 1;
        Eigen::Matrix4d transform;
        localTransform(triangle.block<3, 1>(0, 0), triangle.block<3, 1>(0, 1), triangle.block<3, 1>(0, 2), transform);
        Eigen::Matrix4d localTriangle = transform * triangle_4;

        double l3 = localTriangle(0, 1);
        double u3 = localTriangle(0, 2);
        double v3 = localTriangle(1, 2);

        Eigen::Vector4d tempX;
        tempX.block<3, 1>(0, 0) = x;
        tempX[3] = 1;
        Eigen::Vector4d localX = transform * tempX;
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
        Eigen::Vector3d sminus, splus;
        sminus[0] = -((l3-u3)*(l3-u0)+v3*v0)/l1;
        sminus[1] = -(u3*(u3-u0)+v3*(v3-v0))/l2;
        sminus[2] = -u0;
        splus[0] = ((u3-l3)*(u3-u0)+v3*(v3-v0))/l1;
        splus[1] = (u3*u0+v3*v0)/l2;
        splus[2] = l3-u0;

        // eq (4)
        Eigen::Vector3d t0;
        t0[0] = ((u3-l3)*v0+v3*(l3-u0))/l1;
        t0[1] = (v3*u0-u3*v0)/l2;
        t0[2] = v0;

        // eq (5)
        Eigen::Vector3d tplus, tminus;
        tplus[0] = sqrt((u3-u0)*(u3-u0) + (v3-v0)*(v3-v0));
        tplus[1] = sqrt(u0*u0 + v0*v0);
        tplus[2] = sqrt((l3-u0)*(l3-u0) + v0*v0);
        tminus[0] = tplus[2];
        tminus[1] = tplus[0];
        tminus[2] = tplus[1];

        // line 1, pp. 1450
        Eigen::Vector3d R0;
        for(int i = 0; i < 3; i++)
            R0[i] = sqrt(t0[i]*t0[i] + w0*w0);

        //line 2, pp. 1450
        Eigen::Vector3d Rplus, Rminus;
        for(int i = 0; i < 3; i++)
        {
            Rplus[i] = sqrt(tplus[i]*tplus[i] + w0*w0);
            Rminus[i] = sqrt(tminus[i]*tminus[i] + w0*w0);
        }

        // eq (11)
        Eigen::Vector3d f2;
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
        Eigen::Vector3d beta;
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
        double I1 = 0;
        for(int i = 0; i < 3; i++)
            I1 += t0[i]*f2[i];
        I1 -= fabs(w0) * betaSum;
        return I1;
    }
}
