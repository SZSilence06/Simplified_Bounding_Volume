#ifndef MESH_H
#define MESH_H

#include <vector>
#include <jtflib/mesh/io.h>
#include <hjlib/half_edge/half_edge.h>
#include <hjlib/half_edge/container.h>
#include <hjlib/half_edge/operation.h>
#include <hjlib/half_edge/builder.h>

using namespace zjucad::matrix;
using namespace hj::half_edge;
using std::vector;

struct internal_property
{
  struct vert_t {
    vert_t():id_(-1){}
    vert_t(size_t id): id_(id){}

    void assign(const vert_t& x)
    {
      id_ = x.id_; pos_ = x.pos_;
    }

    size_t id_;
    vector<double> pos_;
  };
  struct edge_t {
    edge_t():id_(-1){}
    edge_t(size_t id): id_(id){}

    void assign(const edge_t& x)
    {
      id_ = x.id_; len_ = x.len_;
    }

    size_t id_;
    double len_;
  };
  struct face_t {
    face_t():id_(-1){}

    face_t(size_t id): id_(id){}

    void assign(const face_t& x)
    {
      id_ = x.id_; normal_ = x.normal_;
    }

    size_t id_;
    vector<double> normal_;
  };
};

typedef half_edge_mesh_t<std_vector,internal_property> HalfEdgeMesh;
typedef HalfEdgeMesh::edges_t::const_iterator itrc_edge_t;
typedef HalfEdgeMesh::faces_t::const_iterator itrc_face_t;
typedef HalfEdgeMesh::verts_t::const_iterator itrc_vert_t;
typedef HalfEdgeMesh::verts_t::iterator       itr_vert_t;
typedef HalfEdgeMesh::faces_t::iterator       itr_face_t;
typedef HalfEdgeMesh::edges_t::iterator       itr_edge_t;


class Mesh
{
private:
    const matrixs_t& triangles;
    const matrixr_t& vertices;

    HalfEdgeMesh halfEdgeMesh;

    double averageEdgeLength;

    matrixr_t vn;
    matrixr_t fn;

    std::vector< std::vector<size_t> > neighbourList;

private:
    void buildHalfEdge();

    bool loadOBJ(std::string filePath);

    void computeAverageEdgeLength();

    void computeNormals();

    void buildNeighbourInfo();

public:
    //Mesh();
    Mesh(const matrixr_t& vertices, const matrixs_t& triangles);

    //bool load(std::string filePath);

    void writeOBJ(std::string filePath);

    inline const matrixr_t& getVertices() const { return vertices; }
    inline const matrixs_t& getTriangles() const { return triangles; }
    inline const HalfEdgeMesh& getHalfEdgeMesh() const { return halfEdgeMesh; }
    inline const matrixr_t& getVerticesNormals() const { return vn; }
    inline const matrixr_t& getFaceNormals() const { return fn; }

    inline double getAverageEdgeLength() const { return averageEdgeLength; }

    std::vector<size_t> findOneRingFaces(int vertex) const;
    const std::vector<size_t>& findOneRingVertices(int vertex) const;
    std::vector<size_t> findOneRingVertices(HalfEdgeMesh::verts_t::const_iterator vert) const;

    std::vector<size_t> findXXXRingVertices(int vertex, int xxx) const;

    //Find out the 'vertex' is in how many rings of 'center'.
    int howManyRing(int vertex, int center) const;

    real_t computeFaceArea(size_t face) const;
    real_t computeVolume() const;
    real_t computeVertexAngle(int vertex) const;
};

#endif // MESH_H
