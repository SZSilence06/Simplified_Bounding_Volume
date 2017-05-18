#include "geometry.h"
#include "math.h"
#include <zjucad/matrix/colon.h>
#include <vector>
#include <limits>
#include "External/eigen3.3/Eigen/Dense"
#include <iostream>

namespace WKYLIB {
    //get 2D angle from 2d coordinate. Range from [0, 2*PI).
    double get_2d_angle(double x, double y)
    {
        if(x == 0)
        {
            if(y > 0)
            {
                return PI / 2;
            }
            else if(y < 0)
            {
                return 3 * PI / 2;
            }
            return 0;
        }
        if(x > 0)
        {
            if( y >= 0)
            {
                return atan( y / x);
            }
            return atan(y / x) + PI * 2;
        }
        return atan(y / x) + PI;
    }

    //Transform 3D triangle to 2D. A will be the original point.
    void tri3D_to_tri2D(const matrixr_t &a,const matrixr_t &b,const matrixr_t &c,
                        matrixr_t &b1,matrixr_t &c1)
    {
        real_t ab = norm(b-a);
        real_t ac = norm(c-a);
        real_t bc = norm(c-b);
        real_t theta = acos((ab*ab+ac*ac-bc*bc)/(2*ab*ac));

        b1 = matrixr_t(2,1);
        c1 = matrixr_t(2,1);

        b1[0] = ab;
        b1[1] = 0;
        c1[0] = ac * cos(theta);
        c1[1] = ac * sin(theta);
    }

    //Transform 2D barycenter to coordinates
    void bary2D_to_coor2D(const matrixr_t &a,const matrixr_t &b,const matrixr_t &c,
                          const matrixr_t &bary,matrixr_t &coor)
    {
        coor = bary[0]*a + bary[1]*b + bary[2]*c;
    }

    //Test whether AB and CD have intersections.
    //if ray == true, then cd is a ray whose direction is from c to d
    //return 0: no intersection;
    //return 1: intersect in the lines;
    //return 2: intersect on point A
    //return 3: intersect on point B
    //return 4: intersect on point C
    //return 5: intersect on point D
    int intersect(const matrixr_t &a,const matrixr_t &b,
                   const matrixr_t &c,const matrixr_t &d,
                   bool ray)
    {
        matrixr_t ab = b - a;
        matrixr_t ac = c - a;
        matrixr_t cd = d - c;
        matrixr_t temp(2,2);
        temp(colon(),0) = ab;
        temp(colon(),1) = -cd;

        real_t det = temp(0,0) * temp(1,1) - temp(0,1) * temp(1,0);
        if(ZERO(det)){
            return 0;
        }
        matrixr_t inv(2,2);
        matrixr_t T;
        inv(0,0) = temp(1,1);
        inv(0,1) = -temp(0,1);
        inv(1.0) = -temp(1,0);
        inv(1,1) = temp(0,0);
        inv /= det;
        T =  inv * ac;
        if(ray){
            if(T[0]<0 || T[0] > 1 || T[1]<0){
                return 0;
            }
        }
        else{
            if(T[0]<0 || T[0] > 1 || T[1]<0 || T[1]>1){
                return 0;
            }
        }
        if(ZERO(T[0])){
            return 2;
        }
        if(ZERO(T[0]-1)){
            return 3;
        }
        if(ZERO(T[1])){
            return 4;
        }
        if(!ray && ZERO(T[1]-1)){
            return 5;
        }
        return 1;
    }

    //Test whether a and b are on same side of the given curve.
    bool is_on_same_side(const matrixr_t &curve,const matrixr_t &a,const matrixr_t &b){
        int cninter = 0,ret,i;
        if((ret = intersect(a,b,curve(colon(),1),curve(colon(),0),true)) == true){
            cninter++;
            if(ret == 4){          //on the end point
                if(curve.size() == 2){
                    return false;
                }
                cninter--;
            }
        }
        for(i=1;i<curve.size(2)-2;i++){
            if((ret = intersect(a,b,curve(colon(),i),curve(colon(),i+1),false)) == true){
                cninter++;
                if(ret == 5){
                    cninter--;
                }
            }
        }
        if(curve.size(2) >2 && ((ret = intersect(a,b,curve(colon(),curve.size(2)-2),curve(colon(),curve.size(2)-1),true)) == true)){
            cninter++;
        }
        if(cninter %2 ){
            return false;
        }
        return true;
    }

    //Test whether point p is on seg AB
    //return 0: not on seg;
    //return 1: on point A;
    //return 2: on point B;
    //return 3: inside the seg;
    int is_on_seg(const matrixr_t &a,const matrixr_t &b,const matrixr_t p)
    {
        matrixr_t ab = b - a;
        matrixr_t ap = p - a;
        matrixr_t bp = p - b;
        if(norm(ap) == 0){
            return 1;
        }
        if(norm(bp) == 0){
            return 2;
        }
        real_t t = ap[0] / ab[0];
        for(int i = 1; i < ap.size(); i++)
        {
            double temp = ap[i] / ab[i];
            if(EQUAL(t, temp) == false)
            {
                return 0;
            }
        }
        if(ZERO(t)){
            return 1;
        }
        if(ZERO(t-1)){
            return 2;
        }
        if(t<0 || t>1){
            return 0;
        }
        return 3;
    }

    //Test whether point p is inside a poly.
    //Note: If p is on the edge or vertex, this function also returns true
    bool is_inside_poly(const matrixr_t &poly,const matrixr_t &p){
        matrixr_t a(2,1);
        a[0] = 999999;
        a[1] = p[1];

        int cnInter = 0;
        for(int i=0;i<poly.size(2)-1;i++){
            if(is_on_seg(poly(colon(),i),poly(colon(),i+1),p)){
                return true;
            }
            int inter;
            if((inter = intersect(p,a,poly(colon(),i),poly(colon(),i+1),false)) != 0){
                cnInter++;
                if(inter == 5){       //intersects on the tail of the edge
                    cnInter--;
                }
            }
        }

        if(is_on_seg(poly(colon(),poly.size(2)-1),poly(colon(),0),p)){
            return true;
        }
        int inter;
        if((inter = intersect(p,a,poly(colon(),poly.size(2)-1),poly(colon(),0),false)) != 0){
            cnInter++;
            if(inter == 5){       //intersects on the tail of the edge
                cnInter--;
            }
        }

        if(cnInter % 2){
            return true;
        }
        return false;
    }

    //Tool function used by is_inside_tetra();
    bool same_side(const matrixr_t &v1, const matrixr_t &v2, const matrixr_t &v3,
                                 const matrixr_t& v4, const matrixr_t& p) {
        matrixr_t normal = cross(v2 - v1, v3 - v1);
        double dotV4 = dot(normal, v4 - v1);
        double dotP = dot(normal, p - v1);
        if(fabs(dotV4) < 1e-3 )
            dotV4 = 0;
        if(fabs(dotP) < 1e-3)
            dotP = 0;
        return (dotV4 >0 && dotP >0) || (dotV4 <0 && dotP <0);
    }

    //Test whether a point is inside a tetrahedron.
    bool is_inside_tetra(const matrixr_t& point, const matrixr_t& tetra)
    {
        return same_side(tetra(colon(), 0), tetra(colon(), 1), tetra(colon(), 2), tetra(colon(), 3), point) &&
               same_side(tetra(colon(), 1), tetra(colon(), 2), tetra(colon(), 3), tetra(colon(), 0), point) &&
               same_side(tetra(colon(), 2), tetra(colon(), 3), tetra(colon(), 0), tetra(colon(), 1), point) &&
               same_side(tetra(colon(), 3), tetra(colon(), 0), tetra(colon(), 1), tetra(colon(), 2), point);
    }

    //Compute area of a triangle
    real_t compute_area(const matrixr_t &a, const matrixr_t &b, const matrixr_t &c)
    {
        matrixr_t ab = b - a;
        matrixr_t ac = c - a;
        return 0.5 * (norm(cross(ab,ac)));
    }

    //Tool function used by compute_volume();
    double SignedVolumeOfTriangle(const matrixr_t &a, const matrixr_t &b, const matrixr_t &c) {
        real_t v321 = c[0] * b[1] * a[2];
        real_t v231 = b[0] * c[1] * a[2];
        real_t v312 = c[0] * a[1] * b[2];
        real_t v132 = a[0] * c[1] * b[2];
        real_t v213 = b[0] * a[1] * c[2];
        real_t v123 = a[0] * b[1] * c[2];
        return (1.0 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123);
    }

    //Compute volume of a triangle mesh
    real_t compute_volume(const matrixr_t &V, const matrixs_t &T)
    {
        real_t sum = 0;
        for(int i = 0; i < T.size(2); i++)
        {
            sum += SignedVolumeOfTriangle(V(colon(),T(0,i)), V(colon(),T(1,i)), V(colon(),T(2,i)));
        }
        return fabs(sum);
    }

    //Compute barycenter coordinates of the point p on trinangle abc.
    //return 1 if p is inside the triangle, and 0 instead.
    int barycentric(const Eigen::Vector3d &a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                    const Eigen::Vector3d &p, Eigen::Vector3d &bary)
    {
        const Eigen::Vector3d u = b - a;
        const Eigen::Vector3d v = c - a;
        const Eigen::Vector3d w = p - a;
        const Eigen::Vector3d vw = v.cross(w);
        const Eigen::Vector3d vu = v.cross(u);
        const Eigen::Vector3d uw = u.cross(w);
        const double denom = 1.0 / vu.norm();
        bary[1] = vu.dot(vw) >= 0 ? vw.norm() * denom : -vw.norm() * denom;
        bary[2] = uw.dot(vu) <= 0 ? uw.norm() * denom : -uw.norm() * denom;
        bary[0] = 1 - bary[1] - bary[2];

        bool result = (bary[0] > 0 || fabs(bary[0]) < 1e-6)
                && (bary[1] >0 || fabs(bary[1]) < 1e-6)
                && (bary[2] > 0 || fabs(bary[2]) < 1e-6);
        return result;
    }

    //Compute barycenter coordinates of the point p on 2d trinangle.
    //return 1 if p is inside the triangle, and 0 instead.
    int barycentric_2D(const Eigen::Vector2d &a, const Eigen::Vector2d& b, const Eigen::Vector2d& c,
                       const Eigen::Vector2d &p, Eigen::Vector3d &bary)
    {
        Eigen::Vector3d a_3d, b_3d, c_3d, p_3d;
        a_3d[0] = a[0];
        a_3d[1] = a[1];
        a_3d[2] = 0;
        b_3d[0] = b[0];
        b_3d[1] = b[1];
        b_3d[2] = 0;
        c_3d[0] = c[0];
        c_3d[1] = c[1];
        c_3d[2] = 0;
        p_3d[0] = p[0];
        p_3d[1] = p[1];
        p_3d[2] = 0;

        return barycentric(a_3d, b_3d, c_3d, p_3d, bary);
    }

    //Tool function used by barycentric_tetra()
    double ScTP(const matrixr_t &a, const matrixr_t &b, const matrixr_t &c)
    {
        // computes scalar triple product
        return dot(a, cross(b, c));
    }

    //Compute barycenter coordinates of the point p on tetrahedron.
    //return 1 if p is inside the tetrahedron, and 0 instead.
    int barycentric_tetra(const matrixr_t &point, const matrixr_t &tetra, matrixr_t &bary)
    {
        matrixr_t vap = point - tetra(colon(), 0);
        matrixr_t vbp = point - tetra(colon(), 1);

        matrixr_t vab = tetra(colon(), 1) - tetra(colon(), 0);
        matrixr_t vac = tetra(colon(), 2) - tetra(colon(), 0);
        matrixr_t vad = tetra(colon(), 3) - tetra(colon(), 0);

        matrixr_t vbc = tetra(colon(), 2) - tetra(colon(), 1);
        matrixr_t vbd = tetra(colon(), 3) - tetra(colon(), 1);
        // ScTP computes the scalar triple product
        double va6 = ScTP(vbp, vbd, vbc);
        double vb6 = ScTP(vap, vac, vad);
        double vc6 = ScTP(vap, vad, vab);
        double vd6 = ScTP(vap, vab, vac);
        double v6 = 1 / ScTP(vab, vac, vad);

        bary = matrixr_t(4, 1);
        bary[0] = va6 * v6;
        bary[1] = vb6 * v6;
        bary[2] = vc6 * v6;
        bary[3] = vd6 * v6;

        if(fabs(bary[0] + bary[1] + bary[2] + bary[3] - 1) < 1e-3 && bary[0] > 0 && bary[1] > 0 && bary[2] > 0)
        {
            return 1;
        }
        return 0;
    }

    //compute intersection between a ray and a triangle abc.
    //return the distance. If no intersection, return -1.
    real_t compute_intersection_triangle_ray(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c,
                                          const Eigen::Vector3d &p, const Eigen::Vector3d &dir)
    {
        Eigen::Vector3d ab = b - a;
        Eigen::Vector3d ac = c - a;
        Eigen::Vector3d plane_normal = ab.cross(ac);
        plane_normal = plane_normal / plane_normal.norm();

        Eigen::Vector3d ap = p - a;
        real_t cos_theta = (ap / ap.norm()).dot(plane_normal);
        real_t dist_p_to_plane = ap.norm() * cos_theta;
        Eigen::Vector3d proj_p = p - dist_p_to_plane * plane_normal;

        Eigen::Vector3d pprojp = proj_p - p;
        real_t cos_theta_projp_p_dir = (pprojp / pprojp.norm()).dot(dir / dir.norm());
        real_t dist = fabs(dist_p_to_plane) / cos_theta_projp_p_dir;
        if(dist < 0)
        {
            return -1;
        }
        Eigen::Vector3d intersect_point = p + dist * dir / dir.norm();
        Eigen::Vector3d bary;
        if(barycentric(a, b, c, intersect_point, bary))
        {
            return dist;
        }
        return -1;
    }

    //Generate a ray direction, given the coordinate x,y in screen space.
    //fov is represented in degree.
    //The origin point is in the left bottom corner.
    void generate_ray_dir(const Eigen::Vector3d &eye, const Eigen::Vector3d& lookAt, const Eigen::Vector3d &up,
                          double fov, double x, double y, double screen_width, double screen_height,
                          Eigen::Vector3d &output_dir)
    {
        Eigen::Vector3d front = (lookAt - eye).normalized();
        Eigen::Vector3d up_new = up.normalized();
        Eigen::Vector3d right = front. cross(up_new);

        double sx = x / screen_width;
        double sy = y /screen_height;
        double fovScale = tan(fov * 0.5 * PI / 180) * 2;

        Eigen::Vector3d r = right * ((sx - 0.5) * fovScale) * screen_width / screen_height;
        Eigen::Vector3d u = up_new * ((sy - 0.5) * fovScale);

        output_dir = r + u + front;
        output_dir = output_dir.normalized();
    }

    //Pick vertex from the screen. If failed to pick, return -1.
    //Param: 'error' shows the maximum distance from the vertex to be picked to the intersection between the
    //         eye ray and the mesh. If the distance is bigger than 'error', the vertex will not be picked.
    int pick_mesh_vertex(const std::vector<Eigen::Vector3d> &vertices, const std::vector<Eigen::Vector3i> &triangles,
                         const Eigen::Vector3d &eye, const Eigen::Vector3d& lookAt, const Eigen::Vector3d &up,
                         double fov, double x, double y, double screen_width, double screen_height,
                         double error)
    {
        Eigen::Vector3d ray_dir;

        generate_ray_dir(eye, lookAt, up, fov, x, y, screen_width, screen_height, ray_dir);

        double min_dist = std::numeric_limits<double>::max();
        for(size_t i = 0; i < triangles.size(); i++)
        {
            double dist = compute_intersection_triangle_ray(vertices[triangles[i][0]],
                                              vertices[triangles[i][1]],
                                              vertices[triangles[i][2]],
                                              eye, ray_dir);
            if(dist < 0 || dist >= min_dist)
            {
                continue;
            }
            min_dist = dist;
        }

        if(min_dist == std::numeric_limits<double>::max())
        {
            return -1;   //The eye ray does not intersect with the mesh.
        }

        Eigen::Vector3d intersection = eye + min_dist * ray_dir;
        int picked_vertex = -1;
        min_dist = std::numeric_limits<double>::max();
        for(size_t i = 0; i < vertices.size(); i++)
        {
            auto& vertex = vertices[i];
            double dist = (vertex - intersection).norm();
            if(dist < min_dist && dist < error)
            {
                min_dist = dist;
                picked_vertex = static_cast<int>(i);
            }
        }
        return picked_vertex;
    }

    //transform barycenter to coordinate
    void bary_to_coor(const matrixr_t& vertices, const matrixs_t& triangles, const matrixr_t& bary,
                      matrixr_t& coor)
    {
        matrixr_t triangle = vertices(colon(), triangles(colon(), static_cast<long int>(bary[2])));

        coor = bary[0] * triangle(colon(), 0) + bary[1] * triangle(colon(), 1)
                + (1 - bary[0] - bary[1]) * triangle(colon(), 2);
    }

    //print a matrix
    void print_matrix(const matrixr_t& a, const std::string& name)
    {
        std::cout << name << " : " << a.size(1) << "*" << a.size(2) << std::endl;
        for(int i = 0; i < a.size(1); i++)
        {
            for(int j = 0; j < a.size(2); j++)
            {
                std::cout << a(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    //print a matrix
    void print_matrix(const matrixs_t& a, const std::string& name)
    {
        std::cout << name << " : " << a.size(1) << "*" << a.size(2) << std::endl;
        for(int i = 0; i < a.size(1); i++)
        {
            for(int j = 0; j < a.size(2); j++)
            {
                std::cout << a(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }
}
