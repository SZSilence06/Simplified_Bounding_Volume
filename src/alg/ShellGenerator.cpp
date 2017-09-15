#include "ShellGenerator.h"
#include "Sampler.h"
#include "Shell.h"
#include "MyPoisson.h"
#include <pcl/io/vtk_io.h>
#include <jtflib/mesh/io.h>
#include <zjucad/matrix/io.h>

using namespace zjucad::matrix;

namespace SBV
{
    ShellGenerator::ShellGenerator(const matrixr_t &vertices, const matrixs_t &triangles)
        : mVertices(vertices),
          mTriangles(triangles)
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

        generateSamples(shell);
        computeDerivative();
        generateOuterShell(shell);

        shell.buildKdTree();
    }

    void ShellGenerator::generateOuterShell(Shell& shell)
    {
        shell.mOuterShell.resize(3, shell.mInnerShell.size(2));
        for(int i = 0; i < shell.mInnerShell.size(2); i++)
        {
            const vec3_t x = shell.mInnerShell(colon(), i);
            shell.mOuterShell(colon(), i) = trace(x);
        }
    }

    vec3_t ShellGenerator::trace(const vec3_t x)
    {
        const double STEP = 0.01;
        double t = 0;
        vec3_t result = x;
        while(t < mDistance)
        {
            vec3_t grad = getGradient(x);
            grad /= norm(grad);
            result -= STEP * grad;
            t += STEP;
        }
        return result;
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

    void ShellGenerator::generateSamples(Shell& shell)
    {
        matrixr_t normals;
        Sampler::poissonDisk(mVertices, mTriangles, mSampleRadius, shell.mInnerShell, normals);

        for(int i = 0; i < shell.mInnerShell.size(2); i++)
        {
            SamplePoint sample;
            sample.position = shell.mInnerShell(colon(), i);
            sample.normal = normals(colon(), i);
            sample.value = 1;
            mSamples.push_back(sample);
            sample.normal = -sample.normal;
            sample.value = -1;
            mSamples.push_back(sample);
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
                else if(isOpposite(mSamples[i], mSamples[j]))
                    B[i] += 0.5 * mSamples[j].value;
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
        for(int i = 0; i < N; i++)
        {
            mSamples[i].derivative = un[i];
        }
    }
}
