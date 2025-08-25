#ifndef UBUFFERS_H
#define UBUFFERS_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/vec3.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <glm/gtc/type_ptr.hpp>

class UniformBuffer {
    private:
    GLuint ubo = 0;
    size_t  size = 0; 
    int bindingPoint = 0;

    public:
    UniformBuffer( size_t  size, int bindingPoint, void* data ) 
    : size( size ), bindingPoint( bindingPoint ) {

        glGenBuffers(1, &ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, ubo);
        glBufferData(GL_UNIFORM_BUFFER, size, nullptr, GL_DYNAMIC_DRAW);
        glBindBufferBase(GL_UNIFORM_BUFFER, bindingPoint, ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);

        update( 0, size, data );
    }

    void update( size_t  offset, size_t  dSize, void *data ) {

        if (offset + dSize > size) {
            std::cerr << "UniformBuffer update out of range!\n";
            return;
        }

        glBindBuffer(GL_UNIFORM_BUFFER, ubo);
        glBufferSubData(GL_UNIFORM_BUFFER, offset, dSize, data);
        glBindBuffer(GL_UNIFORM_BUFFER, 0 );
    }

    void bind() const {
        glBindBufferBase(GL_UNIFORM_BUFFER, bindingPoint, ubo);
    }

    ~UniformBuffer() {
        if (ubo != 0) glDeleteBuffers(1, &ubo);
    }

    GLuint getBindingPoint() const { return bindingPoint; }
    GLuint getBufferID() const { return ubo; }

};

#endif