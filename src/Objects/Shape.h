#ifndef SHAPE_H
#define SHAPE_H

#include <vector>
#include <glm/vec3.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>

typedef struct {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> colors;
    std::vector<glm::vec3> normals;
  } ShapeMesh;  

class Shape {

    private: 
        // std::vector<glm::vec3> vertices, colors, normals;
        ShapeMesh mesh;
        GLuint VAO, VBO;
    public:
        GLenum mode;
        Shape( GLenum mode ): mode( mode ) {
            // glGenBuffers(1, &EBO);
        }
        
        void beginShape() {            
            glGenVertexArrays(1, &VAO);
            glGenBuffers(1, &VBO);
        }

        void endShape() {
            size_t vertices_size = mesh.vertices.size() * sizeof(glm::vec3);
            size_t colors_size  = mesh.colors.size()  * sizeof(glm::vec3);
            size_t normals_size  = mesh.normals.size()  * sizeof(glm::vec3);
            size_t total_size = vertices_size + colors_size + normals_size;

            glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, total_size, nullptr, GL_STATIC_DRAW);
        
            glBufferSubData(GL_ARRAY_BUFFER, 0, vertices_size, mesh.vertices.data());
            glBufferSubData(GL_ARRAY_BUFFER, vertices_size, colors_size, mesh.colors.data());
            glBufferSubData(GL_ARRAY_BUFFER, vertices_size + colors_size, normals_size, mesh.normals.data());

            // vertices at 0 to vertices_size
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
            glEnableVertexAttribArray(0);

            // colors start at vertices_size on
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)vertices_size);
            glEnableVertexAttribArray(1);

            glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)(vertices_size + colors_size));
            glEnableVertexAttribArray(2);
        
            // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            // glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.v_indices.size() * sizeof(unsigned int), mesh.v_indices.data(), GL_STATIC_DRAW);        
        }

        void draw() {
            glBindVertexArray(VAO);
            glDrawArrays(mode, 0, mesh.vertices.size());
            glBindVertexArray(0);
        }

        void color( glm::vec3 col ) {
            mesh.colors.push_back( col );
        }

        void vertex( glm::vec3 v ) {
            mesh.vertices.push_back( v );
        }

        void normal( glm::vec3 n ) {
            mesh.normals.push_back( n );
        }
};

#endif 