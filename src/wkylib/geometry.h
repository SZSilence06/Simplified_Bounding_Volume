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
    const float PI = 3.1415926;

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
    real_t compute_intersection_triangle_ray(const matrixr_t &a, const matrixr_t &b, const matrixr_t &c,
                                          const matrixr_t &p, const matrixr_t &dir);

    //Compute barycenter coordinates of the point p on trinangle.
    //return 1 if p is inside the triangle, and 0 instead.
    //int barycentric(const matrixr_t &point, const matrixr_t &triangle, matrixr_t &bary);

    //Compute barycenter coordinates of the point p on 2d trinangle abc.
    //return 1 if p is inside the triangle, and 0 instead.
    int barycentric_2D(const matrixr_t &point, const matrixr_t &triangle, matrixr_t &bary);

    //Compute barycenter coordinates of the point p on tetrahedron.
    //return 1 if p is inside the tetrahedron, and 0 instead.
    int barycentric_tetra(const matrixr_t &point, const matrixr_t &tetra, matrixr_t &bary);

    //Generate a ray direction, given the coordinate x,y in screen space.
    void generate_ray_dir(const matrixr_t &eye, const matrixr_t& lookAt, const matrixr_t &up,
                          double fov, double x, double y, double screen_width, double screen_height,
                          matrixr_t &output_dir);

    //Pick vertex from the screen. If failed to pick, return -1.
    //Param: 'error' shows the maximum distance from the vertex to be picked to the intersection between the  \
             eye ray and the mesh. If the distance is bigger than 'error', the vertex will not be picked.
    int pick_mesh_vertex(const matrixr_t &vertices, const matrixs_t &triangles, const matrixr_t &eye,
                         const matrixr_t& lookAt, const matrixr_t &up, double fov,
                         double x, double y, double screen_width, double screen_height,
                         double error);

    //transform barycenter to coordinate
    void bary_to_coor(const matrixr_t& vertices, const matrixs_t& triangles, const matrixr_t& bary,
                      matrixr_t& coor);

    template <typename MAT>
    int barycentric(const MAT &point, const MAT &triangle, matrixr_t &bary)
    {
        if(point.size(1) != 3 || triangle.size(1) != 3)
        {
            throw std::invalid_argument("You can only use 3d matrices in barycentric().");
        }

        /*const matrixr_t& a = triangle(colon(), 0);
        const matrixr_t& b = triangle(colon(), 1);
        const matrixr_t& c = triangle(colon(), 2);

        const matrixr_t u = b - a;
        const matrixr_t v = c - a;
        const matrixr_t w = point - a;
        const matrixr_t vw = cross(v, w);
        const matrixr_t vu = cross(v, u);
        const matrixr_t uw = cross(u, w);
        const double denom = 1.0 / norm(vu);
        bary[1] = dot(vu, vw) >= 0 ? norm(vw) * denom : -norm(vw) * denom;
        bary[2] = dot(uw, vu) <= 0 ? norm(uw) * denom : -norm(uw) * denom;
        bary[0] = 1 - bary[1] - bary[2];*/

        Eigen::Vector3d a, b, c, p;
        a[0] = triangle(0, 0);
        a[1] = triangle(1, 0);
        a[2] = triangle(2, 0);
        b[0] = triangle(0, 1);
        b[1] = triangle(1, 1);
        b[2] = triangle(2, 1);
        c[0] = triangle(0, 2);
        c[1] = triangle(1, 2);
        c[2] = triangle(2, 2);
        p[0] = point[0];
        p[1] = point[1];
        p[2] = point[2];

        const Eigen::Vector3d u = b - a;
        const Eigen::Vector3d v = c - a;
        const Eigen::Vector3d w = p - a;
        const Eigen::Vector3d vw = v.cross(w);
        const Eigen::Vector3d vu = v.cross(u);
        const Eigen::Vector3d uw = u.cross(w);
        bary.resize(3, 1);
        double vuNorm = vu.norm();
        if(vuNorm < 1e-8)
            return 0;
        const double denom = 1.0 / vuNorm;
        bary[1] = vu.dot(vw) >= 0 ? vw.norm() * denom : -vw.norm() * denom;
        bary[2] = uw.dot(vu) <= 0 ? uw.norm() * denom : -uw.norm() * denom;
        bary[0] = 1 - bary[1] - bary[2];

        bool result = (bary[0] > 0 || fabs(bary[0]) < 1e-6)
                && (bary[1] >0 || fabs(bary[1]) < 1e-6)
                && (bary[2] > 0 || fabs(bary[2]) < 1e-6);
        return result;
    }
}

#endif // GEOMETRY_H

