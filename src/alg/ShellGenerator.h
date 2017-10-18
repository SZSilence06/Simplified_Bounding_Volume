#ifndef WKY_SHELL_GENERATOR_H
#define WKY_SHELL_GENERATOR_H

#include "Common.h"

namespace pcl
{
    struct PCLPointCloud2;
    struct PolygonMesh;
}

namespace SBV
{
    class Shell;

    class ShellGenerator
    {
    public:
        ShellGenerator(const matrixr_t& vertices, const matrixs_t& triangles, const std::string& outputDirectory);

        void generate(double distance, double sampleRadius, Shell& shell);

    private:
        //void makePCLCloud(const matrixr_t& points, const matrixr_t& normals, pcl::PCLPointCloud2& cloud);
        //void writeHeader(const matrixr_t &points, pcl::PCLPointCloud2& cloud);
        //void makeZJUMesh(const pcl::PolygonMesh& mesh, matrixr_t& vertices, matrixs_t& triangles);

        //for computing Green Function
        struct SamplePoint
        {
            vec3_t position;
            vec3_t normal;
            double value = 0;
            double derivative = 0;
            double size = 0;   //indicating the size of the triangle which the sample point lies in.
            mat3x3_t tri;      //indicating the triangle which the sample point lies in.
            matrixr_t transform;  //matrix for transforming to local
            matrixr_t invTransform;
        };

        void generateSamples(Shell& shell, matrixr_t& normals);
        void addBoundary(Shell& shell);
        void generateOuterShell(Shell& shell, const matrixr_t& inner_shell_normals);
        vec3_t trace(const vec3_t& x, const vec3_t& n);
        void computeDerivative();
        void buildAABB(const Shell& shell, double& xmax, double& xmin, double& ymax, double& ymin, double& zmax, double& zmin);
        void visualizeField(const Shell& shell);
        bool isOpposite(const SamplePoint& a, const SamplePoint& b);

        double G(const vec3_t& x, const vec3_t& x2) { return 1 / (-4 * PI * distance(x, x2)); }

        double Gn(const vec3_t& x, const vec3_t& x2, const vec3_t& n)
        {
            return dot((x - x2) / (4 * PI * pow(distance(x, x2), 3)), n);
        }

        double kernel(const vec3_t& x, const SamplePoint& sample);

        double getFieldValue(const vec3_t& x)
        {
            double result = 0;
            for(int i = 0; i < mSamples.size(); i++)
            {
                //result += kernel(x, mSamples[i]) * mSamples[i].size;
                //if(i == 859)
                result += kernel(x, mSamples[i]);
            }
            //if(result > 1) result = 1;
            //if(result < 0) result = 0;
            return result;
        }

        vec3_t getGradient(const vec3_t& x);

        double distance(const vec3_t& x, const vec3_t& x2);

        void viewTransform(const vec3_t& eye, vec3_t ux, vec3_t uz, matrixr_t& output);
        void localTransform(const vec3_t& a, const vec3_t& b, const vec3_t& c, matrixr_t& output);
        double integrateOverTriangle(const vec3_t& x, const SamplePoint& point, double& I1, vec3_t& Igrad);

    private:
        const matrixr_t& mVertices;
        const matrixs_t& mTriangles;
        double mDistance;
        double mSampleRadius;
        std::string mOutputDirectory;
        const double PI = acos(-1);

        std::vector<SamplePoint> mSamples;
    };
}

#endif
