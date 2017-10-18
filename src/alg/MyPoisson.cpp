/*#include "MyPoisson.h"
#include <pcl/surface/3rdparty/poisson4/octree_poisson.h>
#include <pcl/surface/3rdparty/poisson4/sparse_matrix.h>
#include <pcl/surface/3rdparty/poisson4/function_data.h>
#include <pcl/surface/3rdparty/poisson4/ppolynomial.h>
#include <pcl/surface/3rdparty/poisson4/multi_grid_octree_data.h>
#include <pcl/surface/3rdparty/poisson4/geometry.h>

void MyPoisson::performReconstruction(pcl::PolygonMesh &output)
{
    pcl::poisson::CoredVectorMeshData mesh;
    pcl::poisson::Point3D<float> center;
    float scale = 1.0f;
    switch (this->getDegree())
    {
    case 1:
    {
      myExecute<1> (mesh, center, scale);
      break;
    }
    case 2:
    {
      myExecute<2> (mesh, center, scale);
      break;
    }
    case 3:
    {
      myExecute<3> (mesh, center, scale);
      break;
    }
    case 4:
    {
      myExecute<4> (mesh, center, scale);
      break;
    }
    case 5:
    {
      myExecute<5> (mesh, center, scale);
      break;
    }
    default:
    {
      PCL_ERROR (stderr, "Degree %d not supported\n", this->getDegree());
    }
    }

    // Write output PolygonMesh
    pcl::PointCloud<pcl::PointXYZ> cloud;
    cloud.points.resize (int (mesh.outOfCorePointCount () + mesh.inCorePoints.size ()));
    pcl::poisson::Point3D<float> p;
    for (int i = 0; i < int (mesh.inCorePoints.size ()); i++)
    {
        p = mesh.inCorePoints[i];
        cloud.points[i].x = p.coords[0]*scale+center.coords[0];
        cloud.points[i].y = p.coords[1]*scale+center.coords[1];
        cloud.points[i].z = p.coords[2]*scale+center.coords[2];
    }

    for (int i = int (mesh.inCorePoints.size ()); i < int (mesh.outOfCorePointCount () + mesh.inCorePoints.size ()); i++)
    {
      mesh.nextOutOfCorePoint (p);
      cloud.points[i].x = p.coords[0]*scale+center.coords[0];
      cloud.points[i].y = p.coords[1]*scale+center.coords[1];
      cloud.points[i].z = p.coords[2]*scale+center.coords[2];
    }
    pcl::toPCLPointCloud2 (cloud, output.cloud);
    output.polygons.resize (mesh.polygonCount ());

    // Write faces
    std::vector<pcl::poisson::CoredVertexIndex> polygon;
    for (int p_i = 0; p_i < mesh.polygonCount (); p_i++)
    {
      pcl::Vertices v;
      mesh.nextPolygon (polygon);
      v.vertices.resize (polygon.size ());

      for (int i = 0; i < static_cast<int> (polygon.size ()); ++i)
        if (polygon[i].inCore )
          v.vertices[i] = polygon[i].idx;
        else
          v.vertices[i] = polygon[i].idx + int (mesh.inCorePoints.size ());

      output.polygons[p_i] = v;
    }
}

template<int Degree>
void MyPoisson::myExecute(pcl::poisson::CoredVectorMeshData &mesh, pcl::poisson::Point3D<float> &center, float &scale)
{
    pcl::poisson::Real iso_value = 0;
    pcl::poisson::TreeNodeData::UseIndex = 1;
    pcl::poisson::Octree<Degree> tree;

    /// TODO OPENMP stuff
    //    tree.threads = Threads.value;
    center.coords[0] = center.coords[1] = center.coords[2] = 0;

    if (this->getSolverDivide() < this->getMinDepth())
    {
      PCL_WARN ("[pcl::Poisson] solver_divide_ must be at least as large as min_depth_: %d >= %d\n", this->getSolverDivide(), this->getMinDepth());
      this->setSolverDivide(this->getMinDepth());
    }
    if (this->getIsoDivide()< this->getMinDepth())
    {
      PCL_WARN ("[pcl::Poisson] iso_divide_ must be at least as large as min_depth_: %d >= %d\n", this->getIsoDivide(), this->getMinDepth());
      this->setIsoDivide(this->getMinDepth());
    }

    pcl::poisson::TreeOctNode::SetAllocator (MEMORY_ALLOCATOR_BLOCK_SIZE);

    int kernel_depth_ = this->getDepth() - 2;

    tree.setBSplineData (this->getDepth(), pcl::poisson::Real (1.0 / (1 << this->getDepth())), true);

    tree.maxMemoryUsage = 0;


    int point_count = tree.setTree (this->getInputCloud(), this->getDepth(), this->getMinDepth(), kernel_depth_, this->getSamplesPerNode(),
                                    this->getScale(), center, scale, this->getConfidence(), this->getPointWeight(), true);

    tree.ClipTree ();
    tree.finalize ();
    tree.RefineBoundary (this->getIsoDivide());

    PCL_DEBUG ("Input Points: %d\n" , point_count );
    PCL_DEBUG ("Leaves/Nodes: %d/%d\n" , tree.tree.leaves() , tree.tree.nodes() );

    tree.maxMemoryUsage = 0;
    tree.SetLaplacianConstraints ();

    tree.maxMemoryUsage = 0;
    const int min_iterations = 8;
    const float solver_accuracy = 1e-3f;
    tree.LaplacianMatrixIteration (this->getSolverDivide(), false, min_iterations, solver_accuracy);

    iso_value = tree.GetIsoValue ();
    iso_value += iso_delta_;

    tree.GetMCIsoTriangles (iso_value, this->getIsoDivide(), &mesh, 0, 1, this->getManifold(), this->getOutputPolygons());
}*/
