#ifndef SHADERS_H
#define SHADERS_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/vec3.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <glm/gtc/type_ptr.hpp>

class Shaders {
    private:
    GLuint fs, vs, shader;

    inline void checkCompile(GLuint shader) {
        int success;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            char log[512];
            glGetShaderInfoLog(shader, 512, nullptr, log);
            std::cerr << "Shader compile error: " << log << std::endl;
        }
    }

    public:
    Shaders( const char **vertexShaderSrc, const char **fragmentShaderSrc ) {
        vs = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vs, 1, vertexShaderSrc, nullptr);
        glCompileShader(vs);
        checkCompile(vs);
    
        fs = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fs, 1, fragmentShaderSrc, nullptr);
        glCompileShader(fs);
        checkCompile(fs);
    
        shader = glCreateProgram();
        glAttachShader(shader, vs);
        glAttachShader(shader, fs);
        glLinkProgram(shader);
    }

    GLuint GetShaderID() { return shader; }
};

#endif