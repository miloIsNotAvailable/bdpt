#ifndef SCENE_H
#define SCENE_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/vec3.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <glm/gtc/type_ptr.hpp>

const char* vertexShaderSrc = R"(
    #version 330 core
    layout (location = 0) in vec3 aPos;
    layout (location = 1) in vec3 aColor;
    layout (location = 2) in vec3 aNormal;
    
    out vec4 pos;
    out vec3 vColor;
    out vec3 vNormal;
    
    uniform mat4 uProj;
    uniform mat4 uView;
    
    void main() {
    
        vColor = aColor;
        vNormal = aNormal;

        pos = uProj * uView * vec4(aPos, 1.0);
        gl_Position = pos;
    }
    )";
    
const char* fragmentShaderSrc = R"(
    #version 330 core
    
    in vec4 pos;
    in vec3 vColor;
    in vec3 vNormal;
    // in vec2 gl_PointCoord;

    uniform vec2 uResolution;
    
    uniform sampler2D uTexture;
    
    out vec4 FragColor;
    
    void main() {

        float dist = length(gl_PointCoord - vec2(0.5));
        if (length(gl_PointCoord) > 0. && dist > 0.5)
            discard;

        // FragColor = texture(uTexture, vUV);
        FragColor = vec4( vColor, 1. );
    }
)";


inline void updateEyeFromAngles(glm::vec3& eye, const glm::vec3& center, float radius, float yaw, float pitch) {
    eye.x = center.x + radius * cosf(pitch) * sinf(yaw);
    eye.y = center.y + radius * sinf(pitch);
    eye.z = center.z + radius * cosf(pitch) * cosf(yaw);
}

inline void checkCompile(GLuint shader) {
    int success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char log[512];
        glGetShaderInfoLog(shader, 512, nullptr, log);
        std::cerr << "Shader compile error: " << log << std::endl;
    }
}

class Scene
{
private:
    glm::mat4 perspectiveMatrix;
    glm::mat4 viewMatrix;
    GLuint fs, vs, shader;
    GLint projLoc, viewLoc, res;
    double lastMouseX, lastMouseY;
    
    void updateViewMatrix() {
        updateEyeFromAngles(eye, center, radius, yaw, pitch);
        viewMatrix = glm::lookAt(eye, center, up);
    }

public:
    float w, h;
    float yaw, pitch, fov, aspect, near, far, radius;
    glm::vec3 eye, center, up;
    Scene(float w, float h) : w(w), h(h) {

        fov = 45.0f * 3.14159265f / 180.0f;
        aspect = (float)w / (float)h;
        near = 0.1f;
        far = 100.0f;
        perspectiveMatrix = glm::perspective(fov, aspect, near, far);
        
        eye = glm::vec3(0.0f, 0.0f, 3.0f);        
        center = glm::vec3(0.0f, 0.0f, 0.0f);     
        up = glm::vec3(0.0f, 1.0f, 0.0f);        
        viewMatrix = glm::lookAt(eye, center, up);

        radius = sqrtf(
            (eye.x - center.x)*(eye.x - center.x) + 
            (eye.y - center.y)*(eye.y - center.y) + 
            (eye.z - center.z)*(eye.z - center.z));
        
        yaw = atan2f(eye.x - center.x, eye.z - center.z);
        pitch = asinf((eye.y - center.y) / radius);
    
        updateEyeFromAngles(eye, center, 2., yaw, pitch);

        vs = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vs, 1, &vertexShaderSrc, nullptr);
        glCompileShader(vs);
        checkCompile(vs);
    
        fs = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fs, 1, &fragmentShaderSrc, nullptr);
        glCompileShader(fs);
        checkCompile(fs);
    
        shader = glCreateProgram();
        glAttachShader(shader, vs);
        glAttachShader(shader, fs);
        glLinkProgram(shader);

        res = glGetUniformLocation(shader, "uResolution");
        viewLoc = glGetUniformLocation(shader, "uView");
        projLoc  = glGetUniformLocation(shader, "uProj");

        glDeleteShader(vs);
        glDeleteShader(fs);
    }

    void orbitControls( GLFWwindow * win ) {
        double mouseX, mouseY;
        static float t = 0.;
        static float delta = 0.;
        static float startRadius = radius;
        glfwGetCursorPos(win, &mouseX, &mouseY);
        if( glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS ) {
            // double currentMouseX=mouseX; 
            // double currentMouseY=mouseY;
            glfwGetCursorPos(win, &mouseX, &mouseY);
        
            float deltaX = float(mouseX - lastMouseX);
            float deltaY = float(mouseY - lastMouseY);
        
            const float sensitivity = 0.005f;  // tune this for speed
        
            yaw   -= deltaX * sensitivity;
            pitch += deltaY * sensitivity;
        
            // clamp pitch to avoid flipping upside down
            if (pitch >  1.5f) pitch =  1.5f;  // ~85 degrees
            if (pitch < -1.5f) pitch = -1.5f;
        
            updateViewMatrix();       
        }
        lastMouseX = mouseX;
        lastMouseY = mouseY;

        glfwSetWindowUserPointer(win, this);
        glfwSetScrollCallback(win, [](GLFWwindow* window, double xoffset, double yoffset) {
            auto* self = static_cast<Scene*>(glfwGetWindowUserPointer(window));
            if (self) {
                // self->scroll( window, yoffset );
                delta = yoffset * 1.5;
                startRadius = self->radius;
                t = 0.;
            } else {
                std::cout << "shit done fucked up no user pointer" << std::endl;
            }
        });

        if( t < 1. ) {

            t += 0.1;
            // if( t> 1. ) t = 1.;
        }
        
        float eased = ease( 0., delta, t );
        radius = startRadius - eased;

        updateViewMatrix();
    }   

    // float radOld = radius;
    void scroll( GLFWwindow * win, double yoff ) {
        // radOld = radius;
        radius -= ease( 0, (yoff) * .5, 0.3 );
        updateViewMatrix();
    }

    float ease( float a, float b, float t ) {
        t = sqrt(t);
        return a + (b - a) * t;
    }

    void draw() {
        glClearColor(0.1f, 0.1f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shader);
        // glActiveTexture(GL_TEXTURE0);
        // glUniform1i(glGetUniformLocation(shader, "uTexture"), 0);
        
        glUniform2f(res, w, h);
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(perspectiveMatrix));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(viewMatrix));
    }

    ~Scene() {
        glDeleteProgram(shader);
    };
};

#endif