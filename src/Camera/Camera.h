#ifndef CAMERA_HPP
#define CAMERA_HPP

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/random.hpp>
#include <glm/gtx/norm.hpp>

class Camera {
public:
    glm::vec3 eye;
    glm::vec3 center;
    glm::vec3 ce;         // eye - center
    glm::vec3 nor;        // normalized view direction
    glm::vec3 up;
    glm::vec3 right;
    glm::vec3 focalPlane;
    glm::vec3 c_sensor;

    float apertureSize;
    float focalDist;

    Camera(
        const glm::vec3& center,
        const glm::vec3& eye,
        float apertureSize,
        float focalDist,
        float c_sensor_dist
    ) : center(center), eye(eye), apertureSize(apertureSize), focalDist(focalDist) {

        ce = eye - center;
        nor = glm::normalize(ce);

        c_sensor = center + nor * c_sensor_dist;

        // default up will be set with setCameraBasis()
        focalPlane = center + nor * focalDist;
    }

    void setCameraBasis(const glm::vec3& worldUp = glm::vec3(0.0f, 1.0f, 0.0f)) {
        right = glm::normalize(glm::cross(center - eye, worldUp));
        up = glm::normalize(glm::cross(center - eye, right));
    }

    glm::vec2 sampleRndPointOnDisk(float aperture) {
        float r = glm::linearRand(0.0f, aperture);
        float a = glm::linearRand(0.0f, glm::two_pi<float>());
        return glm::vec2(r * cos(a), r * sin(a));
    }

    glm::vec3 screenToWorldSpace(const glm::vec2& v) {
        return up * v.x + right * v.y;
    }

    float getFocalPointFactor(const glm::vec3& focalPlane, const glm::vec3& primaryRay) {
        glm::vec3 F_c = focalPlane - center;
        float dirDotNor = glm::dot(primaryRay, nor);
        return glm::dot(F_c, nor) / dirDotNor;
    }

    glm::vec3 getFocalPlane(float fDist) const {
        return center + nor * fDist;
    }
};

#endif