#include "Mesh.h"
#include <string.h>
#include <geometry.h>
#include <queue>

namespace WKYLIB {
    Mesh::Mesh(const matrixr_t &vertices, const matrixs_t &triangles)
        : vertices(vertices),triangles(triangles)
    {
        buildHalfEdge();
        computeAverageEdgeLength();
        computeNormals();
        buildNeighbourInfo();
    }

    /*bool Mesh::load(std::string filePath)
    {
        int dotPos = filePath.find_last_of('.');
        if(strcmp(filePath.c_str() + dotPos + 1, "obj") == 0)
        {
            return loadOBJ(filePath);
        }
        else
        {
            return false;
        }
    }*/

    void Mesh::computeNormals()
    {
        vn.resize(3, vertices.size(2));
        fn.resize(3, triangles.size(2));

        matrixr_t vnSum = zeros<real_t>(3, vertices.size(2));

        for(int i = 0; i < triangles.size(2); i++){
            matrixr_t a = vertices(colon(), triangles(0,i));
            matrixr_t b = vertices(colon(), triangles(1,i));
            matrixr_t c = vertices(colon(), triangles(2,i));
            fn(colon(),i) = cross(a-b,a-c);
            for(int j = 0; j < 3; j++){
                vnSum(colon(),triangles(j,i)) += fn(colon(),i);
            }
        }
        for(int i = 0; i < vertices.size(2); i++){
            real_t temp = norm(vnSum(colon(),i));
            vn(colon(),i) = vnSum(colon(),i) / temp;
        }
    }

    void Mesh::computeAverageEdgeLength()
    {
        size_t a = 0, b = 0, cnt = 0;
        double sum = 0;
        for (size_t i = 0; i < triangles.size(2); ++i) {
            for (size_t j = 0; j < 3; ++j) {
                a = j;
                b = (j+1)%3;
                sum += norm(vertices(colon(), triangles(a, i)) - vertices(colon(), triangles(b, i)));
                ++cnt;
            }
        }
        this->averageEdgeLength= sum / cnt;
    }

    void Mesh::writeOBJ(std::string filePath)
    {
        jtf::mesh::save_obj(filePath.c_str(), triangles, vertices);
    }

    void Mesh::buildHalfEdge()
    {
        size_t sizev = vertices.size(2);
        size_t sizet = triangles.size(2);

        halfEdgeMesh = build_topology(triangles);

        size_t i = 0;
        for(itrc_vert_t itv = halfEdgeMesh.verts().begin(); itv != halfEdgeMesh.verts().end(); ++itv,++i){
            internal_property::vert_t vertex;
            vertex.id_ = i;
            for(int j=0;j<3;j++){
                vertex.pos_.push_back(vertices(j,i));
            }
            halfEdgeMesh[itv].assign(vertex);
        }

        i = 0;
        for(itrc_face_t itf = halfEdgeMesh.faces().begin(); itf != halfEdgeMesh.faces().end(); ++itf,++i){
            halfEdgeMesh[itf].id_ = i;
        }
    }

    void Mesh::buildNeighbourInfo()
    {
        for(int i = 0; i < vertices.size(2); i++)
        {
            neighbourList.push_back(std::vector<size_t>());
        }

        auto vert = halfEdgeMesh.verts().begin();
        for(; vert != halfEdgeMesh.verts().end(); ++vert)
        {
            neighbourList[vert->id_] = findOneRingVertices(vert);
        }
    }

    std::vector<size_t> Mesh::findOneRingFaces(int vertex) const
    {
        std::vector<size_t> result;
        auto vert = halfEdgeMesh.verts().begin();
        for(; vert != halfEdgeMesh.verts().end(); ++vert)
        {
            if(vert->id_ == vertex)
            {
                break;
            }
        }

        auto edge = vert->edge();
        do{
            result.push_back(edge->face()->id_);
            edge = edge->next()->oppo();
        }while(edge != vert->edge());

        return result;
    }

    std::vector<size_t> Mesh::findOneRingVertices(HalfEdgeMesh::verts_t::const_iterator vert) const
    {
        std::vector<size_t> result;

        auto edge = vert->edge();
        do{
            result.push_back(edge->oppo()->vert()->id_);
            edge = edge->next()->oppo();
        }while(edge != vert->edge());

        return result;
    }

    const std::vector<size_t>& Mesh::findOneRingVertices(int vertex) const
    {
        return neighbourList[vertex];
    }

    std::vector<size_t> Mesh::findXXXRingVertices(int vertex, int xxx) const
    {
        std::vector<size_t> result;

        auto vert = halfEdgeMesh.verts().begin();
        for(; vert != halfEdgeMesh.verts().end(); ++vert)
        {
            if(vert->id_ == vertex)
            {
                break;
            }
        }
    }

    int Mesh::howManyRing(int vertex, int center) const
    {
        typedef struct{
            int vert;
            int depth = 0;
        }VertInfo;

        std::vector<bool> visited;
        for(int i = 0; i < vertices.size(); i++)
        {
            visited.push_back(false);
        }

        if(vertex == center)
        {
            return 0;
        }

        std::queue<VertInfo> q;
        VertInfo v;
        v.vert = center;
        v.depth = 0;
        q.push(v);
        visited[center] = true;

        //BFS
        while(!q.empty())
        {
            VertInfo& info = q.front();
            q.pop();
            auto& neighbours = neighbourList[info.vert];

            for(int i = 0; i < neighbours.size(); i++)
            {
                if(visited[neighbours[i]])
                {
                    continue;
                }
                if(neighbours[i] == vertex)
                {
                    return info.depth + 1;
                }
                VertInfo newInfo;
                newInfo.vert = neighbours[i];
                newInfo.depth = info.depth + 1;
                visited[neighbours[i]] = true;

                q.push(newInfo);
            }
        }

        //not found
        return -1;
    }

    real_t Mesh::computeFaceArea(size_t face) const
    {
        return compute_area(vertices(colon(),triangles(0,face)), vertices(colon(),triangles(1,face)),
                                 vertices(colon(),triangles(2,face)));
    }

    real_t Mesh::computeVolume() const
    {
        return compute_volume(vertices, triangles);
    }

    real_t Mesh::computeVertexAngle(int vertex) const
    {
        double result = 0;
        auto neighbours = findOneRingVertices(vertex);
        for(int i = 0; i < neighbours.size(); i++)
        {
            int prev = neighbours[(i + neighbours.size() - 1) % neighbours.size()];
            matrixr_t a = vertices(colon(), prev);
            matrixr_t b = vertices(colon(), neighbours[i]);
            matrixr_t c = vertices(colon(), vertex);

            matrixr_t ca = a - c;
            matrixr_t cb = b - c;
            double theta = acos(dot(ca, cb) / (norm(ca) * norm(cb)));
            result += theta;
        }
        return result;
    }
}
