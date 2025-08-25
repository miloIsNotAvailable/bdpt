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
    public:
        Mesh mesh;

        Model( Mesh &mesh ): mesh(mesh) {

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