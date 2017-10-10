#include "ShellGenerator.h"
#include "Sampler.h"
#include "Shell.h"
#include "MyPoisson.h"
#include <pcl/io/vtk_io.h>
#include <jtflib/mesh/io.h>
#include <zjucad/matrix/io.h>
#include <Eigen/IterativeLinearSolvers>
#include <wkylib/mesh/IO.h>
#include <wkylib/geometry.h>
#include "vtk.h"
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkGenericCell.h>

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
        const double scale = 2.5;
        const int res = 100;
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

        //output mesh normals for test
        matrixr_t N(3, mVertices.size(2));
        for(int i = 0; i < mTriangles.size(2); i++)
        {
            const vec3_t a = mVertices(colon(), mTriangles(0, i));
            const vec3_t b = mVertices(colon(), mTriangles(1, i));
            const vec3_t c = mVertices(colon(), mTriangles(2, i));
            vec3_t n = cross(b - a, c - a);
            n /= norm(n);
            N(colon(), mTriangles(0, i)) += n;
            N(colon(), mTriangles(1, i)) += n;
            N(colon(), mTriangles(2, i)) += n;
        }
        for(int i = 0; i < N.size(2); i++)
            N(colon(), i) /= norm(N(colon(), i));
        WKYLIB::Mesh::writeMeshAndNormals((mOutputDirectory + "/mesh.obj"), mVertices, mTriangles, N);

        for(int i = 0; i < mTriangles.size(2); i++)
        {
            const vec3_t a = mVertices(colon(), mTriangles(0, i));
            const vec3_t b = mVertices(colon(), mTriangles(1, i));
            const vec3_t c = mVertices(colon(), mTriangles(2, i));
            vec3_t n = cross(b - a, c - a);
            n /= norm(n);

            SamplePoint sample;
            sample.position = (a + b + c) / 3;
            sample.normal = n;
            sample.value = 1;
            sample.size = WKYLIB::compute_area(a, b, c);
            sample.tri = mVertices(colon(), mTriangles(colon(), i));
            mSamples.push_back(sample);
            sample.normal = -sample.normal;
            sample.value = -1;
            mSamples.push_back(sample);
        }

        addBoundary(shell);
    }

    double ShellGenerator::kernel(const vec3_t &x, const SamplePoint &sample)
    {
        const vec3_t& x2 = sample.position;
        const vec3_t& n = sample.normal;
        double result = sample.value * Gn(x, x2, n) - G(x, x2) * sample.derivative;
        return result;
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
        double xCenter = (xmax + xmin) / 2;
        double yCenter = (ymax + ymin) / 2;
        double zCenter = (zmax + zmin) / 2;
        double radius = std::max(std::max(xCenter - xmin, yCenter - ymin), zCenter - zmin);

        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(xCenter, yCenter, zCenter);
        sphereSource->SetRadius(radius);
        sphereSource->SetThetaResolution(50);
        sphereSource->SetPhiResolution(50);
        sphereSource->Update();
        vtkPolyData* polydata = sphereSource->GetOutput();

        matrixr_t bV(3, polydata->GetNumberOfPoints());
        matrixs_t bT(3, polydata->GetNumberOfCells());

        // Write all of the coordinates of the points in the vtkPolyData to the console.
        for(vtkIdType i = 0; i < polydata->GetNumberOfPoints(); i++)
        {
            vec3_t p;
            polydata->GetPoint(i, &p[0]);
            // This is identical to:
            // polydata->GetPoints()->GetPoint(i,p);
            bV(colon(), i) = p;
        }

        for(vtkIdType i = 0; i < polydata->GetNumberOfCells(); i++)
        {
            vtkCell* cell = polydata->GetCell(i);
            vtkIdList* ids = cell->GetPointIds();
            for(int j = 0; j < 3; j++)
                bT(j, i) = ids->GetId(j);
        }

        matrixr_t normals(3, bV.size(2));
        for(int i = 0; i < bT.size(2); i++)
        {
            const vec3_t a = bV(colon(), bT(0, i));
            const vec3_t b = bV(colon(), bT(1, i));
            const vec3_t c = bV(colon(), bT(2, i));
            vec3_t n = cross(b - a, c - a);
            n /= norm(n);
            normals(colon(), bT(0, i)) += n;
            normals(colon(), bT(1, i)) += n;
            normals(colon(), bT(2, i)) += n;
        }
        for(int i = 0; i < normals.size(2); i++)
            normals(colon(), i) /= norm(normals(colon(), i));

        for(int i = 0; i < bT.size(2); i++)
        {
            const vec3_t a = bV(colon(), bT(0, i));
            const vec3_t b = bV(colon(), bT(1, i));
            const vec3_t c = bV(colon(), bT(2, i));
            vec3_t n = -cross(b - a, c - a);
            n /= norm(n);

            SamplePoint sample;
            sample.position = (a + b + c) / 3;
            sample.normal = n;
            sample.value = 0;
            sample.size = WKYLIB::compute_area(a, b, c);
            sample.tri = bV(colon(), bT(colon(), i));
            mSamples.push_back(sample);
        }

        WKYLIB::Mesh::writeMeshAndNormals((mOutputDirectory + "/boundary.obj"), bV, bT, normals);
    }

    bool ShellGenerator::isOpposite(const SamplePoint &a, const SamplePoint &b)
    {
        return (norm(a.position - b.position) < 1e-6) && (norm(a.normal + b.normal) < 1e-6);
    }

    void ShellGenerator::computeDerivative()
    {
        std::cout << "[INFO] computing derivatives..." << std::endl;

        matrixr_t un;

        const size_t N = mSamples.size();

        matrixr_t A = zeros<double>(N, N);
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                if(j == i || isOpposite(mSamples[i], mSamples[j]))
                {
                    const vec3_t a = mSamples[j].tri(colon(), 0);
                    const vec3_t b = mSamples[j].tri(colon(), 1);
                    const vec3_t c = mSamples[j].tri(colon(), 2);
                    const vec3_t o = (a + b + c) / 3;
                    const vec3_t oab = (o + a + b) / 3;
                    const vec3_t oac = (o + a + c) / 3;
                    const vec3_t obc = (o + b + c) / 3;
                    const double Soab = WKYLIB::compute_area(o, a, b);
                    const double Soac = WKYLIB::compute_area(o, a, c);
                    const double Sobc = WKYLIB::compute_area(o, b, c);
                    A(i, j) = -(G(mSamples[i].position, oab) * Soab + G(mSamples[i].position, oac) * Soac + G(mSamples[i].position, obc) * Sobc);
                }
                else
                    A(i, j) = -G(mSamples[i].position, mSamples[j].position) * mSamples[j].size;
            }
        }

        matrixr_t B = zeros<double>(N, 1);
        for(int i = 0; i < N; i++)
        {
            B[i] += mSamples[i].value;
            for(int j = 0; j < N; j++)
            {
                if(j == i || isOpposite(mSamples[i], mSamples[j]))
                {
                    const vec3_t a = mSamples[j].tri(colon(), 0);
                    const vec3_t b = mSamples[j].tri(colon(), 1);
                    const vec3_t c = mSamples[j].tri(colon(), 2);
                    const vec3_t o = (a + b + c) / 3;
                    const vec3_t oab = (o + a + b) / 3;
                    const vec3_t oac = (o + a + c) / 3;
                    const vec3_t obc = (o + b + c) / 3;
                    const double Soab = WKYLIB::compute_area(o, a, b);
                    const double Soac = WKYLIB::compute_area(o, a, c);
                    const double Sobc = WKYLIB::compute_area(o, b, c);
                    B[i] -= mSamples[j].value * (Gn(mSamples[i].position, oab, mSamples[j].normal) * Soab + Gn(mSamples[i].position, oac, mSamples[j].normal) * Soac
                                                 + Gn(mSamples[i].position, obc, mSamples[j].normal) * Sobc);
                }
                else
                    B[i] -= (mSamples[j].value * Gn(mSamples[i].position, mSamples[j].position, mSamples[j].normal)) * mSamples[j].size;
            }
        }

        Eigen::Map<Eigen::MatrixXd> AA(&A.data()[0], A.size(1), A.size(2));
        Eigen::Map<Eigen::VectorXd> BB(&B.data()[0], B.size(1), B.size(2));

        //std::cout << "[INFO] computing transpose..." << std::endl;
        //BB = AA.transpose() * BB;
        //AA = AA.transpose() * AA;
        std::cout << "[INFO] solving un..." << std::endl;
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

    ShellGenerator::EulerAngle ShellGenerator::eulerAngle(const vec3_t &p0, const vec3_t &pz, const vec3_t &px)
    {
        // translation of local system
        vec3_t uz = pz - p0;
        vec3_t ux = px - p0;

        // unit vectors
        uz /= norm(uz);
        ux /= norm(ux);

        // angle phi from direction of cross product
        EulerAngle result;
        result.phi = atan2(uz[0], -uz[1]);

        // angle theta from dot product of Uz and versor of z-axis
        double eta = -uz[0] * sin(result.phi) + uz[1] * cos(result.phi);
        result.theta = acos(uz[2]);
        if(eta > 0)
            result.theta = -result.theta;

        // angle psi from dot product of Ux and versor of rotated x-axis
        vec3_t uxr;
        uxr[0] = cos(result.phi);
        uxr[1] = sin(result.phi);
        uxr[2] = 0;
        eta = -ux[0] * cos(result.theta) * sin(result.phi) + ux[1] * cos(result.theta) * cos(result.phi) + ux[2] * sin(result.theta);
        double cosarg = dot(ux, uxr);

        // set cosarg between -1 and 1
        if(cosarg > 1)
            cosarg  =1;
        if(cosarg < -1)
            cosarg = -1;

        result.psi = acos(cosarg);
        if(eta < 0)
            result.psi = -result.psi;

        return result;
    }
}
