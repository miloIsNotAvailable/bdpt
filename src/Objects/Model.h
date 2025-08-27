#ifndef MODEL_H
#define MODEL_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/vec3.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <string>
#include "Mesh.h"


struct Hit {
    glm::vec3 hitPoint;
    glm::vec3 hitNormal;
    glm::vec3 color;
    MaterialType mtype;

    Hit( glm::vec3 hitPoint, glm::vec3 hitNormal, glm::vec3 c, MaterialType mtype=MaterialType::Diffuse ) 
    : hitPoint(hitPoint), hitNormal(hitNormal), color( c ), mtype( mtype ) {}
};

class Model {

    private:
      GLuint VAO, VBO, EBO;
      std::vector<Triangle> triangles;
    public:
        Mesh mesh;
        BVH *bvh;

        Model( Mesh &mesh ): mesh(mesh) {

            size_t vert_index = 0;
            for (size_t s = 0; s < mesh.shapes.size(); s++) {
                const auto& shape = mesh.shapes[s];
                for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) {
                    triangles.push_back(Triangle(glm::ivec3(
                        vert_index + 0,
                        vert_index + 1,
                        vert_index + 2
                    )));
                    vert_index += 3;
                }
            }

            std::cout << "building BVH ..." << std::endl;
            bvh = new BVH( triangles, mesh.vertex );
            std::cout << "happily built a BVH ..." << std::endl;

            glGenVertexArrays(1, &VAO);
            glGenBuffers(1, &VBO);
            glGenBuffers(1, &EBO);
            
            size_t vertices_size = mesh.mesh.vertices.size() * sizeof(glm::vec3);
            size_t colors_size  = mesh.mesh.colors.size()  * sizeof(glm::vec3);
            size_t normals_size  = mesh.mesh.normals.size()  * sizeof(glm::vec3);
            size_t total_size = vertices_size + colors_size + normals_size;

            glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, total_size, nullptr, GL_STATIC_DRAW);
        
            glBufferSubData(GL_ARRAY_BUFFER, 0, vertices_size, mesh.mesh.vertices.data());
            glBufferSubData(GL_ARRAY_BUFFER, vertices_size, colors_size, mesh.mesh.colors.data());
            glBufferSubData(GL_ARRAY_BUFFER, vertices_size + colors_size, normals_size, mesh.mesh.normals.data());

            // vertices at 0 to vertices_size
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
            glEnableVertexAttribArray(0);

            // colors start at vertices_size on
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)vertices_size);
            glEnableVertexAttribArray(1);

            glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)(vertices_size + colors_size));
            glEnableVertexAttribArray(2);
        
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.mesh.v_indices.size() * sizeof(unsigned int), mesh.mesh.v_indices.data(), GL_STATIC_DRAW);        
        }

        ~Model() {
            glDeleteVertexArrays(1, &VAO);
            glDeleteBuffers(1, &VBO);
        }

        void draw() {
            glBindVertexArray(VAO);
            glDrawElements(GL_TRIANGLES, mesh.mesh.v_indices.size(), GL_UNSIGNED_INT, 0);
            glBindVertexArray(0);
        }

        bool rayIntersectsAABB( const glm::vec3& origin,
                                const glm::vec3& dir,
                                const glm::vec3& minB,
                                const glm::vec3& maxB,
                                float& tmin,
                                float& tmax )
        {
            tmin = -std::numeric_limits<float>::infinity();
            tmax =  std::numeric_limits<float>::infinity();

            for (int i = 0; i < 3; i++) {
            float invD = 1.0f / dir[i];
            float t0 = (minB[i] - origin[i]) * invD;
            float t1 = (maxB[i] - origin[i]) * invD;
            if (invD < 0.0f) std::swap(t0, t1);

            tmin = t0 > tmin ? t0 : tmin;
            tmax = t1 < tmax ? t1 : tmax;
            if (tmax < tmin) return false;
            }
            return true;
        }

        std::optional<Hit> traverseBVH( const BVH* node,
                                        const std::vector<Vertex>& vertices,
                                        glm::vec3& origin,
                                        glm::vec3& dir,
                                        float *tMin=nullptr, float *tMax=nullptr )
        {
        
            float tmin, tmax;
            if( node == nullptr ) return std::nullopt;
            if( !rayIntersectsAABB( origin, dir, node->minVec, node->maxVec, tmin, tmax ) ) return std::nullopt;

            if(  !node->isLeaf ) {
                std::optional left = traverseBVH(node->left, vertices, origin, dir, tMin, tMax);
                std::optional right = traverseBVH(node->right, vertices, origin, dir, tMin, tMax);
                if (left && right) {
                    float tLeft  = glm::length(left->hitPoint - origin);
                    float tRight = glm::length(right->hitPoint - origin);
                    return (tLeft < tRight) ? left : right;
                }
            
                return left ? left : right;
            }

            std::optional<float> closestHit;
            glm::vec3 hitNormal;
            glm::vec3 hitPoint;
            glm:: vec3 color;
            MaterialType hitType = MaterialType::Diffuse;

            // Find closest intersection
            for (Triangle tri : node->triangles) {
                
                // int k = mesh.mesh.v_indices[tri + 0];
                // int l = mesh.mesh.v_indices[tri + 1];
                // int n = mesh.mesh.v_indices[tri + 2];

                glm::vec3 v0 = vertices[tri.idx.x].position;
                glm::vec3 v1 = vertices[tri.idx.y].position;
                glm::vec3 v2 = vertices[tri.idx.z].position;

                glm::vec3 c0 = vertices[tri.idx.x].color;
                glm::vec3 c1 = vertices[tri.idx.y].color;
                glm::vec3 c2 = vertices[tri.idx.z].color;
                
                Triangle3 tr{v0, v1, v2};
                float u, v;
                auto tHit = ray_intersects_triangle(origin, dir, tr, tMin, tMax, &u, &v);

                if (!tHit.has_value()) continue;

                float t = tHit.value();
                if (!closestHit.has_value() || t < closestHit.value()) {
                    closestHit = t;
                    hitNormal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
                    hitPoint = origin + dir * t;

                    color = u * c0 + v * c1 + (1 - u - v) * c2;
                    
                    hitType = vertices[tri.idx.x].matType;

                    if( vertices[tri.idx.x].light )
                        hitType = MaterialType::Light;
                }
            }
    
            if( !closestHit.has_value() )
                return std::nullopt;

            return Hit( hitPoint, hitNormal, color, hitType );
        }

        std::optional<Hit> isect( glm::vec3& origin,
                                  glm::vec3& dir,
                                  float *tMin=nullptr, float *tMax=nullptr ) {
            return traverseBVH( bvh, mesh.vertex, origin, dir, tMin, tMax );
        }

        std::optional<Hit> isectRay( glm::vec3 &ray_origin, 
                                     glm::vec3 &ray_vector, 
                                     const float *tmin = (const float *)nullptr, 
                                     const float *tmax = (const float *)nullptr ) {

            std::optional<float> closestHit;
            glm::vec3 hitNormal;
            glm::vec3 hitPoint;
            glm:: vec3 color;
            MaterialType hitType = MaterialType::Diffuse;

            // Find closest intersection
            for (int tri = 0; tri < mesh.vertex.size(); tri += 3) {
                
                // int k = mesh.mesh.v_indices[tri + 0];
                // int l = mesh.mesh.v_indices[tri + 1];
                // int n = mesh.mesh.v_indices[tri + 2];
                
                glm::vec3 v0 = mesh.vertex[tri + 0].position;
                glm::vec3 v1 = mesh.vertex[tri + 1].position;
                glm::vec3 v2 = mesh.vertex[tri + 2].position;

                glm::vec3 c0 = mesh.vertex[tri + 0].color;
                glm::vec3 c1 = mesh.vertex[tri + 1].color;
                glm::vec3 c2 = mesh.vertex[tri + 2].color;
                
                Triangle3 tr{v0, v1, v2};
                float u, v;
                auto tHit = ray_intersects_triangle(ray_origin, ray_vector, tr, tmin, tmax, &u, &v);

                if (!tHit.has_value()) continue;

                float t = tHit.value();
                if (!closestHit.has_value() || t < closestHit.value()) {
                    closestHit = t;
                    hitNormal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
                    hitPoint = ray_origin + ray_vector * t;

                    color = u * c0 + v * c1 + (1 - u - v) * c2;
                    
                    hitType = mesh.vertex[tri+0].matType;

                    if( mesh.vertex[tri + 0].light )
                        hitType = MaterialType::Light;
                }
            }

            if( !closestHit.has_value() )
                return std::nullopt;

            return Hit( hitPoint, hitNormal, color, hitType );
        }
};

#endif