#include "ShellGenerator.h"
#include "Sampler.h"
#include "Shell.h"
#include "MyPoisson.h"
//#include <pcl/io/vtk_io.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <jtflib/mesh/io.h>
#include <zjucad/matrix/io.h>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <wkylib/mesh/IO.h>
#include <wkylib/geometry.h>
#include "vtk.h"
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkGenericCell.h>

extern "C"
{
#include "External/fastlap/fastlap.h"
}

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

    const double ONE3  = 0.3333333333333;

    inline double Dot_Product(double* V1, double* V2)
    {
        return V1[0]*V2[0]+V1[1]*V2[1]+V1[2]*V2[2];
    }

    void Cross_Product(double* vector1, double* vector2, double* result_vector)
    {
      result_vector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
      result_vector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
      result_vector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
    }

    double normalize(double vector[3])
    {
      double length;
      int i;

      length = sqrt( vector[0]*vector[0]
            + vector[1]*vector[1]
            + vector[2]*vector[2]);

      for (i=0; i<3; i++) vector[i] = vector[i] / length;

      return length;
    }

    char *delcr(char* str)
    {
      int i, j, k;
      for(k = 0; str[k] != '\0'; k++) if(str[k] == '\n') { str[k] = '\0'; break; }
      for(i = 0; str[i] == ' ' || str[i] == '\t'; i++); /* count leading spaces */
      if(i > 0) {
        for(j = 0; str[j+i] != '\0'; j++) str[j] = str[j+i];
        str[j] = '\0';
      }
      for(k--; str[k] == ' ' || str[k] == '\t'; k--) str[k] = '\0';
      return(str);
    }

    void Dcentroid(int shape, double* pc, double* xcout)
    {
      double corner[4][3], X[3], Y[3], Z[3], vertex1[3], vertex3[3];
      double sum, delta, dl, x1, y1, x2, x3, y3, xc, yc;
      int i, j;
      /* Load the corners. */
      for(i=0; i<4; i++) {
          for(j=0; j<3; j++) {
          corner[i][j] = *(pc++);
          }
      }

      /* Use vertex 0 as the origin and get diags and lengths. */
      for(sum=0, i=0; i<3; i++) {
        X[i] = delta = corner[2][i] - corner[0][i];
        sum += delta * delta;
        vertex1[i] = corner[1][i] - corner[0][i];
        if(shape == QUADRILAT) {
          vertex3[i] = corner[3][i] - corner[0][i];
          Y[i] = corner[1][i] - corner[3][i];
        }
        else if(shape == TRIANGLE) {
          vertex3[i] = corner[2][i] - corner[0][i];
          Y[i] = corner[1][i] - corner[0][i];
        }
        else {
          printf("Dcentroid FE: Shape indicator is neither triangle nor quadrilateral");
          exit(0);
        }
      }
      x2 = sqrt(sum);

      /* Z-axis is normal to two diags. */
      Cross_Product(X, Y, Z);
      normalize(X);
      normalize(Z);

      /* Real Y-axis is normal to X and Z. */
      Cross_Product(Z, X, Y);

      /* Project into the panel axes. */
      y1 = Dot_Product(vertex1, Y);
      y3 = Dot_Product(vertex3, Y);
      x1 = Dot_Product(vertex1, X);
      x3 = Dot_Product(vertex3, X);

      yc = ONE3 * (y1 + y3);
      xc = ONE3 * (x2 + ((x1 * y1 - x3 * y3)/(y1 - y3)));

      *(xcout+0) = corner[0][0] + xc * X[0] + yc * Y[0];
      *(xcout+1) = corner[0][1] + xc * X[1] + yc * Y[1];
      *(xcout+2) = corner[0][2] + xc * X[2] + yc * Y[2];
    }

    void testFastlap()
    {
        const int VERTS = 4, DIMEN = 3;
        const int DIRICHLET = 0, NEUMANN = 1;
        FILE *stream;
        char line[BUFSIZ], title[BUFSIZ];
        int linecnt=0;
        char **chkp, *chk, infile[BUFSIZ], hostname[BUFSIZ];
        double strtod();
        double *exact_sol;
        int size = 1152, nlhs, nrhs, numMom = 4, numLev = 4, i, j;
        int cmderr = FALSE;
        long strtol(), clock;
        char *shapechar;
        double *x, *poten, *dbydnpoten, *xcoll, *xnrm, *lhsvect, *rhsvect;
        int *shape, *type, *dtype, *rhstype, *lhstype, *rhsindex, *lhsindex, job, fljob;
        double error, max_diri=0., ave_diri=0., max_neum=0., ave_neum=0.;
        double cnt_diri = 0, cnt_neum = 0;
          /* Set the tolerance and max iterations for GMRES called by fastlap. */
        double tol = 0.0001;
        int maxit = 32;
        int numit;

        const char* file = "sphere.in";

        if((stream = fopen(file, "r")) == NULL) {
           fprintf(stderr, "\ndriverc FE: Can't open `%s' to read panel data.\n",
               file);
           exit(0);
        }

          time(&clock);
          fprintf(stdout, " Date: %s", ctime(&clock));
          if(gethostname(hostname, BUFSIZ) != -1)
              fprintf(stdout, " Host: %s\n", hostname);
          else fprintf(stdout, " Host: ? (gethostname() failure)\n");


          printf("  Expansion order selected: %i\n",numMom);
          printf("  Depth of tree selected: %i\n",numLev);

          /* Allocate space for the panel vertices and boundary conditions. */
          shape = (int*)calloc(size,sizeof(int));
          x = (double*)calloc(size*VERTS*DIMEN,sizeof(double));
          poten = (double*)calloc(size,sizeof(double));
          dbydnpoten = (double*)calloc(size,sizeof(double));
          type = (int*)calloc(size,sizeof(int));

          /* Allocate space for fastlap arg list vectors. */
          xcoll = (double*)calloc(size*DIMEN,sizeof(double));
          xnrm = (double*)calloc(size*DIMEN,sizeof(double));
          dtype = (int*)calloc(size,sizeof(int));
          lhsvect = (double*)calloc(size,sizeof(double));
          rhsvect = (double*)calloc(size,sizeof(double));
          rhstype = (int*)calloc(size,sizeof(int));
          lhstype = (int*)calloc(size,sizeof(int));
          rhsindex = (int*)calloc(size*VERTS,sizeof(int));
          lhsindex = (int*)calloc(size*VERTS,sizeof(int));

          /* Allocate space for the exact solution vector for error assessment. */
          exact_sol = (double*)calloc(size,sizeof(double));

          /* read in panel data from a file. */
          fgets(line, sizeof(line), stream);
          strcpy(title, delcr(&line[1]));
          while(fgets(line, sizeof(line), stream) != NULL) {
            i = linecnt;
            if(linecnt > (size-1)) {
              fprintf(stderr, "\n More panels than asked for! \n");
              exit(0);
            }
            if(line[0] == 'Q' || line[0] == 'q') {
              shape[i] = QUADRILAT;
              if(sscanf(line,
                        "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
                        &shapechar,
                &x[i*VERTS*DIMEN],&x[i*VERTS*DIMEN+1],&x[i*VERTS*DIMEN+2],
                &x[i*VERTS*DIMEN+3],&x[i*VERTS*DIMEN+4],&x[i*VERTS*DIMEN+5],
                &x[i*VERTS*DIMEN+6],&x[i*VERTS*DIMEN+7],&x[i*VERTS*DIMEN+8],
                &x[i*VERTS*DIMEN+9],&x[i*VERTS*DIMEN+10],&x[i*VERTS*DIMEN+11],
                &poten[i],&dbydnpoten[i],&type[i])
                 != 16) {
                fprintf(stderr, "Bad quad format, line %d:\n%s\n",
                        linecnt, line);
                exit(0);
              }
            }
            else if(line[0] == 'T' || line[0] == 't') {
              shape[i] = TRIANGLE;
              if(sscanf(line,
                        "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
                        &shapechar,
                &x[i*VERTS*DIMEN],&x[i*VERTS*DIMEN+1],&x[i*VERTS*DIMEN+2],
                &x[i*VERTS*DIMEN+3],&x[i*VERTS*DIMEN+4],&x[i*VERTS*DIMEN+5],
                &x[i*VERTS*DIMEN+6],&x[i*VERTS*DIMEN+7],&x[i*VERTS*DIMEN+8],
                &poten[i],&dbydnpoten[i],&type[i])
                 != 13) {
                fprintf(stderr, "Bad tri format, line %d:\n%s\n",
                        linecnt, line);
                exit(0);
              }
            }
            linecnt++;
          }

          printf("  Data file title: %s\n",title);
          printf("  Lines read: %i\n",linecnt);
          size = linecnt;

          /* This is a Green formulation. */
          fljob = 1;
          /* Set up for the fastlap call and save the exact solution for comparison with the
             computed solution.  Note that recovery of the correct signs for Green's Thm. are
             obtained by kidding fastlap about the signs on the lhs and rhs vectors. */
          for(i=0; i<size; i++) {
            if(type[i] == DIRICHLET) {
              rhstype[i] = CONSTANT_DIPOLE;
              lhstype[i] = CONSTANT_SOURCE;
              exact_sol[i] = dbydnpoten[i];
              rhsvect[i] = -poten[i];
              rhsindex[i*VERTS] = i;
              lhsindex[i*VERTS] = i;
              Dcentroid(shape[i], &x[i*VERTS*DIMEN], &xcoll[i*DIMEN]);
              /* fprintf(stdout, "Panel:%d    Centroid:%.8g %.8g %.8g\n",i, xcoll[i*DIMEN],xcoll[i*DIMEN+1],xcoll[i*DIMEN+2]); */
            }
            else if(type[i] == NEUMANN) {
              rhstype[i] = CONSTANT_SOURCE;
              lhstype[i] = CONSTANT_DIPOLE;
              exact_sol[i] = poten[i];
              rhsvect[i] = dbydnpoten[i];
              rhsindex[i*VERTS] = i;
              lhsindex[i*VERTS] = i;
              Dcentroid(shape[i], &x[i*VERTS*DIMEN], &xcoll[i*DIMEN]);
              /* fprintf(stdout, "Panel:%d    Centroid:%.8g %.8g %.8g\n",i, xcoll[i*DIMEN],xcoll[i*DIMEN+1],xcoll[i*DIMEN+2]); */
            }
            else {
              printf("driverc FE: You're missing a boundary condition type");
              exit(0);
            }
          }
          numit = fastlap(&size,&size,&size,x,shape,dtype,lhstype,rhstype,lhsindex,rhsindex,lhsvect,rhsvect,xcoll,xnrm,&numLev,&numMom,&maxit,&tol,&fljob);

          fprintf(stdout, "\n\n %d iterations knocked down residual to:%.8g\n",
              numit, tol);
          /* Compute the average and maximum errors on the Neumann and Dirichlet
             surfaces. Note again, the sign manipulation. */
          double max_error_percent = 0;
          for(i=0;i<size;i++) {
            if(type[i] == DIRICHLET) {
              lhsvect[i] = -lhsvect[i];
              error = sqrt((exact_sol[i] - lhsvect[i])
                   *(exact_sol[i] - lhsvect[i]));
              max_diri = MAX(max_diri,error);
              double error_percent = error / fabs(exact_sol[i]);
              max_error_percent = max_error_percent > error_percent ? max_error_percent : error_percent;
              ave_diri += error;
              cnt_diri += 1.;
              /* fprintf(stdout, "Panel:%d  exact, computed:%.8g %.8g \n",i,exact_sol[i],lhsvect[i]); */
            }
            else if(type[i] == NEUMANN) {
              error = sqrt((exact_sol[i] - lhsvect[i])
                   *(exact_sol[i] - lhsvect[i]));
              max_neum = MAX(max_neum,error);
              ave_neum += error;
              cnt_neum += 1.;
              /* fprintf(stdout, "Panel:%d  exact, computed:%.8g %.8g \n",i,exact_sol[i],lhsvect[i]); */
            }
          }
          if(cnt_diri != 0) {
              ave_diri /= cnt_diri;
              fprintf(stdout, "\nAverage absolute error on Dirichlet surface =%.8g\n",
                  ave_diri);
              fprintf(stdout, "Maximum absolute error on Dirichlet surface =%.8g\n",
                  max_diri);
          }
          if(cnt_neum != 0) {
              ave_neum /= cnt_neum;
              fprintf(stdout, "\nAverage absolute error on Neumann surface =%.8g\n",
                  ave_neum);
              fprintf(stdout, "Maximum absolute error on Neumann surface =%.8g\n",
                  max_neum);
          }

          free(shape);
          free(x);
          free(poten);
          free(dbydnpoten);
          free(type);
          free(xcoll);
          free(xnrm);
          free(dtype);
          free(lhsvect);
          free(rhsvect);
          free(lhstype);
          free(rhstype);
          free(lhsindex);
          free(rhsindex);
          free(exact_sol);

          return;
    }

    void ShellGenerator::generate(double distance, double sampleRadius, Shell &shell)
    {
        mSampleRadius = sampleRadius;
        mDistance = distance;

        matrixr_t normals;

        //testFastlap();
        //exit(0);

        generateSamples(shell, normals);
        //computeDerivative();
        computeDerivative_fastlap();
        visualizeField(shell, true);
        exit(0);
        generateOuterShell(shell, normals);
        //exit(0);

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
    void grid2vtk(OS &os, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int res, bool planar) {
        const int resX = planar ? 1 : res;
        const int resY = res;
        const int resZ = res;

      os << "# vtk DataFile Version 2.0\nSample rectilinear grid\nASCII\nDATASET RECTILINEAR_GRID\n";
      os << "DIMENSIONS" << " " << resX << " " << resY << " " << resZ << std::endl;

      // order: first x, then y, finally z
      const double dx = (xmax - xmin)/(res-1);
      const double dy = (ymax - ymin)/(res-1);
      const double dz = (zmax - zmin)/(res-1);

      os << "X_COORDINATES " << resX << " float\n";
      if(planar)
      {
          for (size_t i = res / 2 - 1; i < res / 2; ++i)
              os << xmin+i*dx << std::endl;
      }
      else
      {
          for (size_t i = 0; i < res; ++i)
              os << xmin+i*dx << std::endl;
      }

      os << "Y_COORDINATES " << res << " float\n";
      for (size_t i = 0; i < res; ++i)
        os << ymin+i*dy << std::endl;

      os << "Z_COORDINATES " << res << " float\n";
      for (size_t i = 0; i < res; ++i)
        os << zmin+i*dz << std::endl;
    }

    void ShellGenerator::visualizeField(const Shell& shell, bool planar)
    {
        std::cout << "[INFO] Generating field..." << std::endl;
        const double scale = 2.5;
        const int res = 200;
        double xmax, xmin, ymax, ymin, zmax, zmin;
        buildAABB(shell, xmax, xmin, ymax, ymin, zmax, zmin);
        scaleAABB(xmin, xmax, ymin, ymax, zmin, zmax, scale);

        std::ofstream of(mOutputDirectory + "/field.vtk");
        grid2vtk(of, xmin, xmax, ymin, ymax, zmin, zmax, res, planar);

        const double xStep = (xmax - xmin) / (res - 1);
        const double yStep = (ymax - ymin) / (res - 1);
        const double zStep = (zmax - zmin) / (res - 1);

        const int resX = planar ? 1 : res;
        const int resY = res;
        const int resZ = res;
        matrixr_t fieldData(resX * resY * resZ, 1);

#pragma omp parallel for
        for(size_t i = 0; i < resZ; i++) {
            for(size_t j = 0; j < resY; j++) {
                if(planar) {
                    for(size_t k = res / 2 - 1; k < res / 2; k++) {
                        const size_t idx = i * resY * resX + j * resX + 0;
                        vec3_t pos;
                        pos[0] = xmin + k * xStep;
                        pos[1] = ymin + j * yStep;
                        pos[2] = zmin + i * zStep;
                        double value = getFieldValue(pos);
                        fieldData[idx] = value;
                    }
                }
                else {
                    for(size_t k = 0; k < resX; k++) {
                        const size_t idx = i * resY * resX + j * resX + k;
                        vec3_t pos;
                        pos[0] = xmin + k * xStep;
                        pos[1] = ymin + j * yStep;
                        pos[2] = zmin + i * zStep;
                        double value = getFieldValue(pos);
                        fieldData[idx] = value;
                    }
                }
            }
        }

        point_data(of, fieldData.begin(), fieldData.size(), "field");
        of.close();
        std::cout << "[INFO] field generated." << std::endl;
    }

    void ShellGenerator::generateOuterShell(Shell& shell, const matrixr_t& inner_shell_normals)
    {
        std::cout << "[INFO] generating outer shell..." << std::endl;
        shell.mOuterShell.resize(3, shell.mInnerShell.size(2));
#pragma omp parallel for
        for(int i = 0; i < shell.mInnerShell.size(2); i++)
        {
            const vec3_t x = shell.mInnerShell(colon(), i);
            shell.mOuterShell(colon(), i) = trace(x, inner_shell_normals(colon(), i));
        }
        std::cout << "[INFO] outer shell generated." << std::endl;
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

    double ShellGenerator::getFieldValue(const vec3_t &x)
    {
        double result = 0;
        for(int i = 0; i < mSamples.size(); i++)
        {
            result += kernel(x, mSamples[i]);
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
        const double STEP = 0.0001;
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

        WKYLIB::Mesh::writePointsAndNormals(mOutputDirectory + "/inner_shell_normal.vtk", shell.mInnerShell, normals);

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

        //for test fastlap
        std::ifstream in("value.dat");
        std::vector<double> values;
        for(int i = 0; i < mTriangles.size(2); i++) {
            double read;
            in >> read;
            values.push_back(read);
        }
        in.close();

        for(int i = 0; i < mTriangles.size(2); i++)
        {
            const vec3_t a = mVertices(colon(), mTriangles(0, i));
            const vec3_t b = mVertices(colon(), mTriangles(1, i));
            const vec3_t c = mVertices(colon(), mTriangles(2, i));
            vec3_t n = cross(b - a, c - a);
            n /= norm(n);

            SamplePoint sample;
            sample.position = (a + b + c) / 3;
            sample.normal = -n;
            sample.value = 1;
            sample.size = WKYLIB::compute_area(a, b, c);
            sample.tri = mVertices(colon(), mTriangles(colon(), i));
            localTransform(sample.tri(colon(), 0), sample.tri(colon(), 1), sample.tri(colon(), 2), sample.transform);
            sample.invTransform = sample.transform;
            inv(sample.invTransform);
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

    void ShellGenerator::addBoundary(const Shell& shell)
    {
        //boundary is a scaled bounding sphere
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
            vec3_t n = cross(b - a, c - a);
            n /= norm(n);

            SamplePoint sample;
            sample.position = (a + b + c) / 3;
            sample.normal = n;
            sample.value = 0;
            sample.size = WKYLIB::compute_area(a, b, c);
            sample.tri = bV(colon(), bT(colon(), i));
            localTransform(sample.tri(colon(), 0), sample.tri(colon(), 1), sample.tri(colon(), 2), sample.transform);
            sample.invTransform = sample.transform;
            inv(sample.invTransform);
            mSamples.push_back(sample);
        }

        WKYLIB::Mesh::writeMeshAndNormals((mOutputDirectory + "/boundary.obj"), bV, bT, normals);

        /*const double INF = 100000;
        matrixr_t bV(3, 8);
        bV(0, 0) = -INF; bV(1, 0) = -INF; bV(2, 0) = -INF;
        bV(0, 1) = INF; bV(1, 1) = -INF; bV(2, 1) = -INF;
        bV(0, 2) = -INF; bV(1, 2) = INF; bV(2, 2) = -INF;
        bV(0, 3) = INF; bV(1, 3) = INF; bV(2, 3) = -INF;
        bV(0, 4) = -INF; bV(1, 4) = -INF; bV(2, 4) = INF;
        bV(0, 5) = INF; bV(1, 5) = -INF; bV(2, 5) = INF;
        bV(0, 6) = -INF; bV(1, 6) = INF; bV(2, 6) = INF;
        bV(0, 7) = INF; bV(1, 7) = INF; bV(2, 7) = INF;

        matrixs_t bT(3, 12);
        bT(0, 0) = 0; bT(1, 0) = 2; bT(2, 0) = 1;
        bT(0, 1) = 1; bT(1, 1) = 2; bT(2, 1) = 3;
        bT(0, 2) = 0; bT(1, 2) = 4; bT(2, 2) = 6;
        bT(0, 3) = 0; bT(1, 3) = 6; bT(2, 3) = 2;
        bT(0, 4) = 0; bT(1, 4) = 1; bT(2, 4) = 5;
        bT(0, 5) = 0; bT(1, 5) = 5; bT(2, 5) = 4;
        bT(0, 6) = 6; bT(1, 6) = 4; bT(2, 6) = 5;
        bT(0, 7) = 6; bT(1, 7) = 5; bT(2, 7) = 7;
        bT(0, 8) = 7; bT(1, 8) = 5; bT(2, 8) = 1;
        bT(0, 9) = 7; bT(1, 9) = 1; bT(2, 9) = 3;
        bT(0, 10) = 2; bT(1, 10) = 6; bT(2, 10) = 7;
        bT(0, 11) = 2; bT(1, 11) = 7; bT(2, 11) = 3;

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
            vec3_t n = cross(b - a, c - a);
            n /= norm(n);

            SamplePoint sample;
            sample.position = (a + b + c) / 3;
            sample.normal = n;
            sample.value = 0;
            sample.size = WKYLIB::compute_area(a, b, c);
            sample.tri = bV(colon(), bT(colon(), i));
            localTransform(sample.tri(colon(), 0), sample.tri(colon(), 1), sample.tri(colon(), 2), sample.transform);
            sample.invTransform = sample.transform;
            inv(sample.invTransform);
            mSamples.push_back(sample);
        }

        WKYLIB::Mesh::writeMeshAndNormals((mOutputDirectory + "/boundary.obj"), bV, bT, normals);*/
    }

    bool ShellGenerator::isOpposite(const SamplePoint &a, const SamplePoint &b)
    {
        return (norm(a.position - b.position) < 1e-6) && (norm(a.normal + b.normal) < 1e-6);
    }

    void ShellGenerator::computeDerivative()
    {
        const size_t N = mSamples.size();
        std::cout << "[INFO] computing derivatives of " << N << " samples... " << std::endl;

        matrixr_t un;   

        matrixr_t A = zeros<double>(N, N);
        matrixr_t B = zeros<double>(N, 1);
        for(int i = 0; i < N; i++)
        {
            B[i] = mSamples[i].value;
            for(int j = 0; j < N; j++)
            {
                double I1;
                vec3_t Igrad;
                integrateOverTriangle(mSamples[i].position, mSamples[j], I1, Igrad);
                A(i, j) = I1;
                //A(i, j) = I1 / (4 * PI);
                //B[i] += (mSamples[j].value * dot(Igrad, mSamples[j].normal)) / (4 * PI);
            }
        }

        Eigen::Map<Eigen::MatrixXd> AA(&A.data()[0], A.size(1), A.size(2));
        Eigen::Map<Eigen::VectorXd> BB(&B.data()[0], B.size(1), B.size(2));

        std::cout << "[INFO] solving un..." << std::endl;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(AA);
        Eigen::VectorXd UNN = solver.solve(BB);
        un.resize(UNN.rows(), 1);
        for(int i = 0; i < UNN.rows(); i++){
            un[i] = UNN[i];
        }
        std::cout << un << std::endl;
        std::cout << "un solved." << std::endl;
        for(int i = 0; i < N; i++)
        {
            mSamples[i].derivative = un[i];
        }
    }

    void ShellGenerator::computeDerivative_fastlap()
    {
        int size = mSamples.size();
        double* x = new double[12*size];
        for(int i = 0; i < mSamples.size(); i++)
            for(int j = 0; j < 3; j++)
                for(int k = 0; k < 3; k++)
                    x[12*i + 3*j + k] = mSamples[i].tri(k, j);

        int* shape = new int[size];
        for(int i = 0; i < size; i++)
            shape[i] = TRIANGLE;

        int* dtype = new int[size];
        for(int i = 0; i < size; i++)
            dtype[i] = 0;

        int* lhstype = new int[size];
        int* rhstype = new int[size];
        for(int i = 0; i < size; i++)
        {
            lhstype[i] = CONSTANT_SOURCE;
            rhstype[i] = CONSTANT_DIPOLE;
        }

        int* lhsindex = new int[size * 4];
        int* rhsindex = new int[size * 4];
        for(int i = 0; i < size; i++)
        {
            lhsindex[i * 4] = i;
            rhsindex[i * 4] = i;
        }

        double* lhsVect = new double[size];
        double* rhsVect = new double[size];
        for(int i = 0; i < size; i++)
            rhsVect[i] = mSamples[i].value;

        double* xcoll = new double[size * 3];
        for(int i = 0; i < size; i++)
        {
            Dcentroid(TRIANGLE, &x[i*12], &xcoll[i*3]);
        }

        double* xnrm = new double[size * 3];

        int numLev = 4, numMom = 4, maxit = 32;
        double tol = 1e-4;
        int fljob = INDIRECT;

        fastlap(&size, &size, &size, x, shape, dtype, lhstype, rhstype, lhsindex, rhsindex, lhsVect, rhsVect, xcoll, xnrm, &numLev, &numMom, &maxit, &tol, &fljob);

        matrixr_t error(size, 1);
        for(int i = 0; i < size; i++)
            error[i] = lhsVect[i] - mSamples[i].derivative;

        std::cout << "error is " << error << std::endl;

        for(int i = 0; i < size; i++)
        {
            mSamples[i].derivative = lhsVect[i];
        }

        delete[] x;
        delete[] shape;
        delete[] dtype;
        delete[] lhstype;
        delete[] rhstype;
        delete[] lhsindex;
        delete[] rhsindex;
        delete[] lhsVect;
        delete[] rhsVect;
        delete[] xcoll;
        delete[] xnrm;
    }

    void ShellGenerator::viewTransform(const vec3_t &eye, vec3_t ux, vec3_t uz, mat4x4_t &output)
    {
        uz /= norm(uz);
        ux /= norm(ux);
        const vec3_t uy = cross(uz, ux);

        output(0, 0) = ux[0];
        output(0, 1) = ux[1];
        output(0, 2) = ux[2];
        output(0, 3) = -dot(eye, ux);
        output(1, 0) = uy[0];
        output(1, 1) = uy[1];
        output(1, 2) = uy[2];
        output(1, 3) = -dot(eye, uy);
        output(2, 0) = uz[0];
        output(2, 1) = uz[1];
        output(2, 2) = uz[2];
        output(2, 3) = -dot(eye, uz);
        output(3, 0) = 0;
        output(3, 1) = 0;
        output(3, 2) = 0;
        output(3, 3) = 1;
    }

    void ShellGenerator::localTransform(const vec3_t &a, const vec3_t &b, const vec3_t &c, mat4x4_t &output)
    {
        vec3_t ux = b - a;
        vec3_t uz = cross(ux, c - a);
        output.resize(4, 4);
        viewTransform(a, ux, uz, output);
    }

    //closed-form calculation according to Graglia 1993.
    double ShellGenerator::integrateOverTriangle(const vec3_t& x, const SamplePoint &point, double& I1, vec3_t& Igrad)
    {
        mat4x4_t triangle;
        triangle(colon(0, 2), colon(0, 2)) = point.tri;
        triangle(3, colon(0, 2)) = ones<double>(1, 3);
        mat4x4_t localTriangle = point.transform * triangle;

        double l3 = localTriangle(0, 1);
        double u3 = localTriangle(0, 2);
        double v3 = localTriangle(1, 2);

        if(l3 < 0)
            throw std::runtime_error("l3 < 0.");

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
        double threshold = 1e-6 * std::min(std::min(l1,l2), l3);
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
                if(temp < 0)
                    std::cerr << "[WARNING] computing log of negative number. i = " << i
                              << " tplus[0] = " << tplus[0]
                              << " tminus[0] = " << tminus[0]
                              << " splus[0] = " << splus[0]
                              << " sminus[0] = " << sminus[0]
                              << ". line " << __LINE__ << std::endl;
            }
            else
            {
                 temp = (Rplus[i]+splus[i]) / (Rminus[i]+sminus[i]);
                 if(temp < 0)
                     std::cerr << "[WARNING] computing log of negative number. i = " << i << ". line " << __LINE__ << std::endl;
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

    double ShellGenerator::kernel(const vec3_t &x, const SamplePoint &sample)
    {
        double I1;
        vec3_t Igrad;
        integrateOverTriangle(x, sample, I1, Igrad);
        //double result = (sample.value * dot(Igrad, sample.normal) - sample.derivative * I1) / (-4 * PI);
        double result = I1 * sample.derivative;
        return result;
    }
}
