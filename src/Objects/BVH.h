#ifndef BVH_H
#define BVH_H

#include <glm/glm.hpp>
#include <glm/gtx/component_wise.hpp> 
#include <vector>
#include <algorithm>
#include <iostream>
// #include "Mesh.h"

enum class MaterialType {
    Light,
    Diffuse,
    Glossy,
    Mirror,
    Dieletric
};

struct Vertex {
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec3 color;
    MaterialType matType;
    bool light=false;
};


class Triangle {
    public:
    // glm::vec3 u, v, w, c;
    // Triangle( glm::vec3 v0, glm::vec3 v1, glm::vec3 v2 ): u(v0), v(v1), w(v2) {
    //     c = (u + v + w) / 3.f;
    // }
    glm::ivec3 idx;
    Triangle( glm::ivec3 idx ): idx(idx) {}

    glm::vec3 centroid( std::vector<Vertex> v ) {
        glm::vec3 centroid = (v[idx.x].position +
            v[idx.y].position +
            v[idx.z].position) / 3.0f;

            return centroid;
    }
};

class BVH {
    public:
    std::vector<Triangle> triangles;
    BVH* left;
    BVH* right;
    glm::vec3 minVec;
    glm::vec3 maxVec;
    bool isLeaf = false;

    BVH( std::vector<Triangle> &triangles, const std::vector<Vertex>& vertices ) : triangles(triangles) {
        minVec = glm::vec3( std::numeric_limits<float>::infinity() );
        maxVec = glm::vec3( -std::numeric_limits<float>::infinity() );

        aabb( vertices );
    
        if (triangles.size() <= 8) {
            isLeaf = true;
            return;
        }
        
        glm::vec3 extents = maxVec - minVec;
        
        int longestAxisIndex = 0;
        if (extents.y > extents.x) longestAxisIndex = 1;
        if (extents.z > extents[longestAxisIndex]) longestAxisIndex = 2;        
        
        sortByAxis( longestAxisIndex, vertices );
        
        int mid = triangles.size() / 2.;
        std::vector<Triangle> leftTris(triangles.begin(), triangles.begin() + mid);
        std::vector<Triangle> rightTris(triangles.begin() + mid, triangles.end());

        left = new BVH(leftTris, vertices);
        right = new BVH(rightTris, vertices);
        
        triangles.clear();   
    }
      
    void aabb(const std::vector<Vertex>& vertices) {
        for (const auto& t : triangles) {
            const glm::vec3& u = vertices[t.idx.x].position;
            const glm::vec3& v = vertices[t.idx.y].position;
            const glm::vec3& w = vertices[t.idx.z].position;

            minVec = glm::min(minVec, glm::min(u, glm::min(v, w)));
            maxVec = glm::max(maxVec, glm::max(u, glm::max(v, w)));
        }
    }

    void sortByAxis(int axis, const std::vector<Vertex>& vertices) {
        std::sort(triangles.begin(), triangles.end(),
            [&](Triangle& a, Triangle& b) {
                glm::vec3 ac = a.centroid( vertices );
                glm::vec3 bc = b.centroid( vertices );
                return ac[axis] < bc[axis];
            });
    }
};

#endif