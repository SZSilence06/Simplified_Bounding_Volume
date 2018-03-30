#ifndef WKY_SHELL_GENERATOR_H
#define WKY_SHELL_GENERATOR_H

#include "Common.h"
#include "Tracer.h"

namespace pcl
{
    struct PCLPointCloud2;
    struct PolygonMesh;
}

namespace SBV
{
    class Shell;
    class FieldComputer;

    /**
     * @brief The ShellGenerator class
     *
     * This class is used to generate the outer shell using Green's function.
     */
    class ShellGenerator
    {
    public:
        ShellGenerator(const matrixr_t& vertices, const matrixs_t& triangles, const std::string& outputDirectory);

        /**
         * @brief generate shells.
         * @param distance : between the outer shell and the given mesh.
         * @param sampleRadius : poisson sample radius to sample the mesh.
         * @param shell : output shells
         * @param isVisualizeField : true to visualize the field instead of generating shells, and false otherwise
         */
        void generate(double distance, double sampleRadius, Shell& shell, bool isVisualizeField);

    private:
        //void makePCLCloud(const matrixr_t& points, const matrixr_t& normals, pcl::PCLPointCloud2& cloud);
        //void writeHeader(const matrixr_t &points, pcl::PCLPointCloud2& cloud);
        //void makeZJUMesh(const pcl::PolygonMesh& mesh, matrixr_t& vertices, matrixs_t& triangles);

        /**
         * @brief generate samples on the given mesh.
         * @param shell : output shells. After this call, only sample points on inner shell will be generated.
         * @param normals : output normals of the sample points.
         */
        void generateSamples(Shell& shell, matrixr_t& normals);

        // Add the infinity bounding sphere boundary of the mesh.
        void addBoundary(const Shell& shell);

        // generate the outer shell.
        void generateOuterShell(Shell& shell, const matrixr_t& inner_shell_normals);

        /**
         * @brief This function computes the derivatives using direct solver. Now it is deprecated, and is replaced by computeDerivative_fastlap().
         */
        void computeDerivative();

        /**
         * @brief This function computes the derivatives using fastlap.
         */
        void computeDerivative_fastlap();

        // compute AABB of the shell.
        void buildAABB(const Shell& shell, double& xmax, double& xmin, double& ymax, double& ymin, double& zmax, double& zmin);

        /**
         * @brief Visualize the potential field.
         * @param shell : the given shell.
         * @param planar : true to visualize only on a plane, and false to visualize a 3d field.
         */
        void visualizeField(const Shell& shell, bool planar);

        // the green function. It is deprecated, since the computation is done on GPU now.
        double G(const vec3_t& x, const vec3_t& x2) { return 1 / (-4 * PI * distance(x, x2)); }

        // normal derivative of Green function. It is deprecated, since the computation is done on GPU now.
        double Gn(const vec3_t& x, const vec3_t& x2, const vec3_t& n)
        {
            return dot((x - x2) / (4 * PI * pow(distance(x, x2), 3)), n);
        }

        // function to be integrated on triangles. It is deprecated, since the computation is done on GPU now.
        double kernel(const vec3_t& x, const Triangle& sample);

        // get the field value on the given point.
        double getFieldValue(const vec3_t& x);

        // get the gradient on the given point.
        vec3_t getGradient(const vec3_t& x);

        // compute distance between two points.
        double distance(const vec3_t& x, const vec3_t& x2);

        /**
         * @brief compute the view transform matrix to transform from world coordinate system to a local coordinate system.
         * @param eye : original point of the local coordinate system.
         * @param ux : x axis of the local coordinate system.
         * @param uz : z axis of the local coordinate system.
         * @param output : output transform matrix.
         */
        void viewTransform(const vec3_t& eye, vec3_t ux, vec3_t uz, mat4x4_t& output);

        /**
         * @brief compute the view transform matrix to transform from world coordinate system to a local coordinate system.
         * @param a : original point of the local coordinate system.
         * @param b : a point on the x axis of the local coordinate system.
         * @param c : a point on the y axis of the local coordinate system.
         * @param output : output transform matrix.
         */
        void localTransform(const vec3_t& a, const vec3_t& b, const vec3_t& c, mat4x4_t& output);

        /**
         * @brief compute the green function integration on a triangle.
         *        This function is deprecated, since the computation is done on GPU now.
         * @param x : observer point.
         * @param point : the triangle to be integrated.
         * @param I1 : output integrate result of the Green function.
         * @param Igrad : output integrate result of the gradient of the Green function.
         */
         void integrateOverTriangle(const vec3_t& x, const Triangle& point, double& I1, vec3_t& Igrad);

    private:
        const matrixr_t& mVertices;
        const matrixs_t& mTriangles;
        double mDistance;
        double mSampleRadius;
        std::string mOutputDirectory;
        const double PI = acos(-1);

        std::vector<Triangle> mSamples;

        FieldComputer mField;
    };
}

#endif
