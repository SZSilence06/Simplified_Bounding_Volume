/*#ifndef WKY_PCL_MY_POISSON_H
#define WKY_PCL_MY_POISSON_H

#include <pcl/surface/poisson.h>

class MyPoisson : public pcl::Poisson<pcl::PointNormal>
{
public:
    virtual void performReconstruction(pcl::PolygonMesh &output);

    template<int Degree> void
    myExecute (pcl::poisson::CoredVectorMeshData &mesh,
               pcl::poisson::Point3D<float> &center,
               float &scale);

    inline void setIsoDelta(double iso_delta) { iso_delta_ = iso_delta; }

private:
    double iso_delta_ = 0;
};

#endif*/
