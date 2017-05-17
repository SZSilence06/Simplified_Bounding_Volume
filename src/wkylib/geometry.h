#ifndef WKY_GEOMETRY_H
#define WKY_GEOMETRY_H

#include <zjucad/matrix/matrix.h>
#include <eigen3/Eigen/Dense>

using namespace zjucad::matrix;

typedef double real_t;
typedef matrix<real_t> matrixr_t;
typedef matrix<size_t> matrixs_t;

namespace WKYLIB {
    #define ZERO(x) (fabs(x)<1e-6)
    #define EQUAL(a,b) (ZERO(a-b))
    const double PI = 3.1415926;

    //get 2D angle from 2d coordinate. Range from [0, 2*PI).
    double get_2d_angle(double x, double y);

    //Transform 3D triangle to 2D. A will be the original point.
    void tri3D_to_tri2D(const matrixr_t &a,const matrixr_t &b,const matrixr_t &c,
                        matrixr_t &b1,matrixr_t &c1);

    //Transform 2D barycenter to coordinates
    void bary2D_to_coor2D(const matrixr_t &a,const matrixr_t &b,const matrixr_t &c,
                          const matrixr_t &bary,matrixr_t &coor);

    //Test whether AB and CD have intersections.
    //if ray == true, then cd is a ray whose direction is from c to d
    //Note: If AB and CD intersect on the vertex, this function returns 2
    //Note: If AB and CD are parallel,this function returns 0
    int intersect(const matrixr_t &a,const matrixr_t &b,
                   const matrixr_t &c,const matrixr_t &d,
                   bool ray);

    //Test whether a and b are on same side of the given curve.
    bool is_on_same_side(const matrixr_t &curve,const matrixr_t &a,const matrixr_t &b);

    //Test whether point p is on seg AB
    int is_on_seg(const matrixr_t &a,const matrixr_t &b,const matrixr_t p);

    //Test whether point p is inside a poly.
    //Note: If p is on the edge or vertex, this function also returns true
    bool is_inside_poly(const matrixr_t &poly,const matrixr_t &p);

    //Test whether a point is inside a tetrahedron.
    bool is_inside_tetra(const matrixr_t& point, const matrixr_t& tetra);

    //Compute area of a triangle;
    real_t compute_area(const matrixr_t &a, const matrixr_t &b, const matrixr_t &c);

    //Compute volume of a triangle mesh
    real_t compute_volume(const matrixr_t &V, const matrixs_t &T);

    //compute intersection between a ray and a triangle abc.
    //return the distance. If no intersection, return -1.
    real_t compute_intersection_triangle_ray(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c,
                                             const Eigen::Vector3d &p, const Eigen::Vector3d &dir);

    //Compute barycenter coordinates of the point p on trinangle.
    //return 1 if p is inside the triangle, and 0 instead.
    int barycentric(const Eigen::Vector3d &a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                    const Eigen::Vector3d &p, Eigen::Vector3d &bary);

    //Compute barycenter coordinates of the point p on 2d trinangle abc.
    //return 1 if p is inside the triangle, and 0 instead.
    int barycentric_2D(const Eigen::Vector2d &a, const Eigen::Vector2d& b, const Eigen::Vector2d& c,
                       const Eigen::Vector2d &p, Eigen::Vector3d &bary);

    //Compute barycenter coordinates of the point p on tetrahedron.
    //return 1 if p is inside the tetrahedron, and 0 instead.
    int barycentric_tetra(const matrixr_t &point, const matrixr_t &tetra, matrixr_t &bary);

    //Generate a ray direction, given the coordinate x,y in screen space.
    void generate_ray_dir(const matrixr_t &eye, const matrixr_t& lookAt, const matrixr_t &up,
                          double fov, double x, double y, double screen_width, double screen_height,
                          matrixr_t &output_dir);

    //Pick vertex from the screen. If failed to pick, return -1.
    //Param: 'error' shows the maximum distance from the vertex to be picked to the intersection between the
    //         eye ray and the mesh. If the distance is bigger than 'error', the vertex will not be picked.
    int pick_mesh_vertex(const matrixr_t &vertices, const matrixs_t &triangles, const matrixr_t &eye,
                         const matrixr_t& lookAt, const matrixr_t &up, double fov,
                         double x, double y, double screen_width, double screen_height,
                         double error);

    //transform barycenter to coordinate
    void bary_to_coor(const matrixr_t& vertices, const matrixs_t& triangles, const matrixr_t& bary,
                      matrixr_t& coor);
}

#endif // GEOMETRY_H

