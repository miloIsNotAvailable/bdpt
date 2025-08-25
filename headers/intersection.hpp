#ifndef INTERSECTION_H
#define INTERSECTION_H

#include <glm/glm.hpp>
#include <glm/gtc/epsilon.hpp>
#include <optional>
#include <limits>

struct Triangle3 {
    glm::vec3 a;
    glm::vec3 b;
    glm::vec3 c;
};

std::optional<float> ray_intersects_triangle(
    glm::vec3& ray_origin,
    glm::vec3& ray_vector,
    Triangle3& triangle,
    const float *tmin=nullptr,
    const float *tmax=nullptr,
    float *u_=nullptr,
    float *v_=nullptr );

    #endif