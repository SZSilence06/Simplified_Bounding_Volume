#include "ShellGenerator.h"
#include "Sampler.h"
#include "Shell.h"
#include "MyPoisson.h"
#include <pcl/io/vtk_io.h>
#include <jtflib/mesh/io.h>
#include <zjucad/matrix/io.h>
#include <Eigen/IterativeLinearSolvers>
#include <wkylib/mesh/IO.h>
#include "vtk.h"

using namespace zjucad::matrix;

namespace SBV
{
    ShellGenerator::ShellGenerator(const matrixr_t &vertices, const matrixs_t &triangles, const std::string& outputDirectory)
        : mVertices(vertices),
          mTriangles(triangles),
          mOutputDirectory(outputDirectory)
    {

    }

    /*void ShellGenerator::generate(double distance, double sampleRadius, Shell& shell)
    {
        matrixr_t normals;
        Sampler::poissonDisk(mVertices, mTriangles, sampleRadius, shell.mInnerShell, normals);

        pcl::PCLPointCloud2 cloud;
        makePCLCloud(shell.mInnerShell, normals, cloud);

        pcl::PointCloud<pcl::PointNormal>::Ptr xyz_cloud(new pcl::PointCloud<pcl::PointNormal>());
        pcl::fromPCLPointCloud2(cloud, *xyz_cloud);

        pcl::PolygonMesh mesh;
        MyPoisson poisson;
        poisson.setDepth(8);
        poisson.setSolverDivide(8);
        poisson.setIsoDivide(8);
        poisson.setPointWeight(4);
        poisson.setInputCloud(xyz_cloud);
        poisson.setIsoDelta(0.5);
        poisson.reconstruct(mesh);

        matrixr_t recon_vertices;
        matrixs_t recon_triangles;
        makeZJUMesh(mesh, recon_vertices, recon_triangles);
        jtf::mesh::save_obj("test2.obj", recon_triangles, recon_vertices);

        matrixr_t normal_outer;
        Sampler::poissonDisk(recon_vertices, recon_triangles, sampleRadius, shell.mOuterShell, normal_outer);

        shell.buildKdTree();
    }

    void ShellGenerator::makePCLCloud(const matrixr_t &points, const matrixr_t &normals, pcl::PCLPointCloud2 &cloud)
    {
        writeHeader(points, cloud);

        //write xyz
        for (int i = 0; i < points.size(2); i++)
        {
            for(int j = 0; j < 3; j++)
            {
                float value = points(j, i);
                memcpy (&cloud.data[i * cloud.point_step + cloud.fields[j].offset],
                        &value,
                        sizeof (float));
            }
        }

        // Get normal_x fields indices
        int normal_x_field = -1;
        for (std::size_t i = 0; i < cloud.fields.size (); ++i)
            if (cloud.fields[i].name == "normal_x")
            {
                normal_x_field = i;
                break;
            }

        //write normals
        for (int i = 0; i < normals.size(2); i++)
        {
            for(int j = 0; j < 3; j++)
            {
                float value = normals(j, i);
                memcpy (&cloud.data[i * cloud.point_step + cloud.fields[normal_x_field + j].offset],
                        &value,
                        sizeof (float));
            }
        }
    }

    void ShellGenerator::writeHeader(const matrixr_t &points, pcl::PCLPointCloud2 &cloud)
    {
        int field_offset = 0;
        for (int i = 0; i < 3; ++i, field_offset += 4)
        {
            cloud.fields.push_back (pcl::PCLPointField ());
            cloud.fields[i].offset   = field_offset;
            cloud.fields[i].datatype = pcl::PCLPointField::FLOAT32;
            cloud.fields[i].count    = 1;
        }
        cloud.fields[0].name = "x";
        cloud.fields[1].name = "y";
        cloud.fields[2].name = "z";

        std::string normals_names[3] = { "normal_x", "normal_y", "normal_z" };
        for (int i = 0; i < 3; ++i, field_offset += 4)
        {
            cloud.fields.push_back (pcl::PCLPointField ());
            pcl::PCLPointField& last = cloud.fields.back ();
            last.name     = normals_names[i];
            last.offset   = field_offset;
            last.datatype = pcl::PCLPointField::FLOAT32;
            last.count    = 1;
        }

        cloud.point_step = field_offset;
        cloud.width      = points.size(2);
        cloud.height     = 1;
        cloud.row_step   = cloud.point_step * cloud.width;
        cloud.is_dense   = true;
        cloud.data.resize (cloud.point_step * points.size(2));
    }

    void ShellGenerator::makeZJUMesh(const pcl::PolygonMesh &mesh, matrixr_t &vertices, matrixs_t& triangles)
    {
        //number of points
        int nr_points  = mesh.cloud.width * mesh.cloud.height;
        unsigned point_size = static_cast<unsigned> (mesh.cloud.data.size () / nr_points);
        //number of faces
        unsigned nr_faces = static_cast<unsigned> (mesh.polygons.size ());

        vertices.resize(3, nr_points);
        triangles.resize(3, nr_faces);

        //write vertex xyz
        for(int i = 0; i < nr_points; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                float value;
                memcpy (&value, &mesh.cloud.data[i * point_size + mesh.cloud.fields[j].offset], sizeof (float));
                vertices(j, i) = value;
            }
        }

        //write triangle faces
        for(int i = 0; i < nr_faces; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                triangles(j, i) = mesh.polygons[i].vertices[j];
            }
        }
    }*/

    void ShellGenerator::generate(double distance, double sampleRadius, Shell &shell)
    {
        mSampleRadius = sampleRadius;
        mDistance = distance;

        matrixr_t normals;

        generateSamples(shell, normals);
        computeDerivative();
        visualizeField(shell);
        generateOuterShell(shell, normals);

        shell.buildKdTree();
    }

    void scaleAABB(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax, double scale)
    {
        double xCenter = (xmin + xmax) / 2;
        double yCenter = (ymin + ymax) / 2;
        double zCenter = (zmin + zmax) / 2;
        double xLength = xmax - xmin;
        double yLength = ymax - ymin;
        double zLength = zmax - zmin;
        xLength *= scale;
        yLength *= scale;
        zLength *= scale;
        xmin = xCenter - xLength / 2;
        xmax = xCenter + xLength / 2;
        ymin = yCenter - yLength / 2;
        ymax = yCenter + yLength / 2;
        zmin = zCenter - zLength / 2;
        zmax = zCenter + zLength / 2;
    }

    template <typename OS>
    void grid2vtk(OS &os, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int res) {
      os << "# vtk DataFile Version 2.0\nSample rectilinear grid\nASCII\nDATASET RECTILINEAR_GRID\n";
      os << "DIMENSIONS" << " " << res << " " << res << " " << res << std::endl;

      // order: first x, then y, finally z
      const double dx = (xmax - xmin)/(res-1);
      const double dy = (ymax - ymin)/(res-1);
      const double dz = (zmax - zmin)/(res-1);

      os << "X_COORDINATES " << res << " float\n";
      for (size_t i = 0; i < res; ++i)
        os << xmin+i*dx << std::endl;

      os << "Y_COORDINATES " << res << " float\n";
      for (size_t i = 0; i < res; ++i)
        os << ymin+i*dy << std::endl;

      os << "Z_COORDINATES " << res << " float\n";
      for (size_t i = 0; i < res; ++i)
        os << zmin+i*dz << std::endl;
    }

    void ShellGenerator::visualizeField(const Shell& shell)
    {
        const double scale = 3;
        const int res = 50;
        double xmax, xmin, ymax, ymin, zmax, zmin;
        buildAABB(shell, xmax, xmin, ymax, ymin, zmax, zmin);
        scaleAABB(xmin, xmax, ymin, ymax, zmin, zmax, scale);

        std::ofstream of(mOutputDirectory + "/field.vtk");
        grid2vtk(of, xmin, xmax, ymin, ymax, zmin, zmax, res);

        const double xStep = (xmax - xmin) / res;
        const double yStep = (ymax - ymin) / res;
        const double zStep = (zmax - zmin) / res;
        matrixr_t fieldData(res * res * res, 1);

#pragma omp parallel for
        for(size_t i = 0; i < res; i++) {
            for(size_t j = 0; j < res; j++) {
                for(size_t k = 0; k < res; k++)
                {
                    const size_t idx = i * res * res + j * res + k;
                    vec3_t pos;
                    pos[0] = xmin + k * xStep;
                    pos[1] = ymin + j * yStep;
                    pos[2] = zmin + i * zStep;
                    double value = getFieldValue(pos);
                    //if(value > 1)
                    //    value = 1;
                    //if(value < 0)
                    //    value = 0;
                    fieldData[idx] = value;
                }
            }
        }

        point_data(of, fieldData.begin(), fieldData.size(), "field");
        of.close();
        std::cout << "[INFO] field generated." << std::endl;
    }

    void ShellGenerator::generateOuterShell(Shell& shell, const matrixr_t& inner_shell_normals)
    {
        shell.mOuterShell.resize(3, shell.mInnerShell.size(2));
        for(int i = 0; i < shell.mInnerShell.size(2); i++)
        {
            const vec3_t x = shell.mInnerShell(colon(), i);
            shell.mOuterShell(colon(), i) = trace(x, inner_shell_normals(colon(), i));
        }
    }

    vec3_t ShellGenerator::trace(const vec3_t& x, const vec3_t& n)
    {
        const double STEP = 0.01;
        double t = 0;
        vec3_t result = x;
        while(t < mDistance)
        {
            if(t == 0)
            {
                result += STEP * n;
            }
            else
            {
                vec3_t grad = getGradient(result);
                double norm_grad = norm(grad);
                if(norm_grad == 0) {
                    std::cerr << "[warning] gradient equals zero" << std::endl;
                }
                grad /= norm_grad;
                result -= STEP * grad;
            }
            t += STEP;
        }
        return result;
    }

    double ShellGenerator::distance(const vec3_t &x, const vec3_t &x2)
    {
        const double r = norm(x - x2);
        if(fabs(r) < 1e-6) {
            std::cerr << "# warning: near the boundary " << r << std::endl;
        }
        return r;
    }

    vec3_t ShellGenerator::getGradient(const vec3_t &x)
    {
        const double STEP = 0.01;
        double a = getFieldValue(x);

        vec3_t xx = x;
        xx[0] += STEP;
        double b = getFieldValue(xx);

        vec3_t xy = x;
        xy[1] += STEP;
        double c = getFieldValue(xy);

        vec3_t xz = x;
        xz[2] += STEP;
        double d = getFieldValue(xz);

        vec3_t result;
        result[0] = (b - a) / STEP;
        result[1] = (c - a) / STEP;
        result[2] = (d - a) / STEP;
        return result;
    }

    void ShellGenerator::generateSamples(Shell& shell, matrixr_t& normals)
    {
        Sampler::poissonDisk(mVertices, mTriangles, mSampleRadius, shell.mInnerShell, normals);

        for(int i = 0; i < shell.mInnerShell.size(2); i++)
        {
            SamplePoint sample;
            sample.position = shell.mInnerShell(colon(), i);
            sample.normal = normals(colon(), i);
            sample.value = 1;
            mSamples.push_back(sample);
            //sample.normal = -sample.normal;
            //sample.value = -1;
            //mSamples.push_back(sample);
        }

        addBoundary(shell);
    }

    void ShellGenerator::buildAABB(const Shell& shell, double &xmax, double &xmin, double &ymax, double &ymin, double &zmax, double &zmin)
    {
        xmax = std::numeric_limits<double>::lowest();
        xmin = std::numeric_limits<double>::max();
        ymax = std::numeric_limits<double>::lowest();
        ymin = std::numeric_limits<double>::max();
        zmax = std::numeric_limits<double>::lowest();
        zmin = std::numeric_limits<double>::max();
        for(int i = 0; i < shell.mInnerShell.size(2); i++)
        {
            xmax = std::max(xmax, shell.mInnerShell(0, i));
            xmin = std::min(xmin, shell.mInnerShell(0, i));
            ymax = std::max(ymax, shell.mInnerShell(1, i));
            ymin = std::min(ymin, shell.mInnerShell(1, i));
            zmax = std::max(zmax, shell.mInnerShell(2, i));
            zmin = std::min(zmin, shell.mInnerShell(2, i));
        }
    }

    void ShellGenerator::addBoundary(Shell& shell)
    {
        const double scale = 2;
        double xmax, xmin, ymax, ymin, zmax, zmin;
        buildAABB(shell, xmax, xmin, ymax, ymin, zmax, zmin);
        scaleAABB(xmin, xmax, ymin, ymax, zmin, zmax, scale);

        matrixr_t bV(3, 8);
        matrixs_t bT(3, 12);
        bV(0, 0) = xmin; bV(1, 0) = ymin; bV(2, 0) = zmin;
        bV(0, 1) = xmin; bV(1, 1) = ymin; bV(2, 1) = zmax;
        bV(0, 2) = xmin; bV(1, 2) = ymax; bV(2, 2) = zmin;
        bV(0, 3) = xmin; bV(1, 3) = ymax; bV(2, 3) = zmax;
        bV(0, 4) = xmax; bV(1, 4) = ymin; bV(2, 4) = zmin;
        bV(0, 5) = xmax; bV(1, 5) = ymin; bV(2, 5) = zmax;
        bV(0, 6) = xmax; bV(1, 6) = ymax; bV(2, 6) = zmin;
        bV(0, 7) = xmax; bV(1, 7) = ymax; bV(2, 7) = zmax;

        bT(0, 0) = 0; bT(1, 0) = 1; bT(2, 0) = 2;
        bT(0, 1) = 1; bT(1, 1) = 2; bT(2, 1) = 3;
        bT(0, 2) = 0; bT(1, 2) = 1; bT(2, 2) = 4;
        bT(0, 3) = 1; bT(1, 3) = 4; bT(2, 3) = 5;
        bT(0, 4) = 0; bT(1, 4) = 2; bT(2, 4) = 4;
        bT(0, 5) = 2; bT(1, 5) = 4; bT(2, 5) = 6;
        bT(0, 6) = 4; bT(1, 6) = 5; bT(2, 6) = 6;
        bT(0, 7) = 5; bT(1, 7) = 6; bT(2, 7) = 7;
        bT(0, 8) = 2; bT(1, 8) = 3; bT(2, 8) = 6;
        bT(0, 9) = 3; bT(1, 9) = 6; bT(2, 9) = 7;
        bT(0, 10) = 1; bT(1, 10) = 3; bT(2, 10) = 5;
        bT(0, 11) = 3; bT(1, 11) = 5; bT(2, 11) = 7;

        matrixr_t samples, normals;
        Sampler::poissonDisk(bV, bT, mSampleRadius, samples, normals);

        for(int i = 0; i < samples.size(2); i++)
        {
            SamplePoint point;
            point.position = samples(colon(), i);
            point.normal = -normals(colon(), i);
            point.value = 0;
            mSamples.push_back(point);
        }
    }

    bool ShellGenerator::isOpposite(const SamplePoint &a, const SamplePoint &b)
    {
        return (norm(a.position - b.position) < 1e-6) && (norm(a.normal + b.normal) < 1e-6);
    }

    void ShellGenerator::computeDerivative()
    {
        matrixr_t un;

        const double l = mSampleRadius;
        const size_t N = mSamples.size();

        matrixr_t A = zeros<double>(N, N);
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                if(j == i || isOpposite(mSamples[i], mSamples[j]))
                    A(i, j) = -l / 2;
                else
                    A(i, j) = -G(mSamples[i].position, mSamples[j].position) * PI * l * l;
            }
        }

        matrixr_t B = zeros<double>(N, 1);
        for(int i = 0; i < N; i++)
        {
            B[i] += mSamples[i].value;
            for(int j = 0; j < N; j++)
            {
                if(j == i)
                    B[i] -= 0.5 * mSamples[j].value;
                //    B[i] -= 0;
                else if(isOpposite(mSamples[i], mSamples[j]))
                    B[i] += 0.5 * mSamples[j].value;
                //    B[i] += 0;
                else
                    B[i] -= (mSamples[j].value * Gn(mSamples[i].position, mSamples[j].position, mSamples[j].normal)) * PI * l * l;
            }
        }

        std::cout << "solving un..." << std::endl;
        Eigen::Map<Eigen::MatrixXd> AA(&A.data()[0], A.size(1), A.size(2));
        Eigen::Map<Eigen::VectorXd> BB(&B.data()[0], B.size(1), B.size(2));

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(AA);
        Eigen::VectorXd UNN = solver.solve(BB);
        un.resize(UNN.rows(), 1);
        for(int i = 0; i < UNN.rows(); i++){
            un[i] = UNN[i];
        }
        std::cout << "un solved." << std::endl;
        std::cout << un << std::endl;
        for(int i = 0; i < N; i++)
        {
            mSamples[i].derivative = un[i];
        }
    }
}
