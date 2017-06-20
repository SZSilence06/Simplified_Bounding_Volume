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
        ShellGenerator(const matrixr_t& vertices, const matrixs_t& triangles);

        void generate(double distance, double sampleRadius, Shell& shell);

    private:
        void makePCLCloud(const matrixr_t& points, const matrixr_t& normals, pcl::PCLPointCloud2& cloud);
        void writeHeader(const matrixr_t &points, pcl::PCLPointCloud2& cloud);
        void makeZJUMesh(const pcl::PolygonMesh& mesh, matrixr_t& vertices, matrixs_t& triangles);

    private:
        const matrixr_t& mVertices;
        const matrixs_t& mTriangles;
    };
}

#endif
