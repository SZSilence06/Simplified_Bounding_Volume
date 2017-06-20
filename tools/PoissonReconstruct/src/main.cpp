#define CGAL_EIGEN3_ENABLED

#include <pcl/PCLPointCloud2.h>
#include <pcl/io/obj_io.h>
#include <pcl/surface/poisson.h>
#include "MyPoisson.h"

int main(int argc, char** argv)
{
    pcl::PCLPointCloud2 cloud;
    pcl::io::loadOBJFile(argv[1], cloud);

    pcl::PolygonMesh mesh;
    pcl::PointCloud<pcl::PointNormal>::Ptr xyz_cloud(new pcl::PointCloud<pcl::PointNormal>());
    pcl::fromPCLPointCloud2(cloud, *xyz_cloud);

    MyPoisson poisson;
    poisson.setDepth(8);
    poisson.setSolverDivide(8);
    poisson.setIsoDivide(8);
    poisson.setPointWeight(4);
    poisson.setInputCloud(xyz_cloud);
    poisson.reconstruct(mesh);

    pcl::io::saveOBJFile("reconstruction.obj", mesh);

    return 0;
}
