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
    float *v_=nullptr )
{
    constexpr float epsilon = std::numeric_limits<float>::epsilon();

    glm::vec3 edge1 = triangle.b - triangle.a;
    glm::vec3 edge2 = triangle.c - triangle.a;
    glm::vec3 ray_cross_e2 = glm::cross(ray_vector, edge2);
    float det = glm::dot(edge1, ray_cross_e2);

    if (det > -epsilon && det < epsilon)
        return std::nullopt; // Parallel ray

    float inv_det = 1.0f / det;
    glm::vec3 s = ray_origin - triangle.a;
    float u = inv_det * glm::dot(s, ray_cross_e2);
    *u_ = u;

    if ((u < 0 && glm::abs(u) > epsilon) || (u > 1 && glm::abs(u - 1) > epsilon))
        return std::nullopt;

    glm::vec3 s_cross_e1 = glm::cross(s, edge1);
    float v = inv_det * glm::dot(ray_vector, s_cross_e1);
    *v_ = v;

    if ((v < 0 && glm::abs(v) > epsilon) || (u + v > 1 && glm::abs(u + v - 1) > epsilon))
        return std::nullopt;

    float t = inv_det * glm::dot(edge2, s_cross_e1);

    if( (tmin && tmax) && (t > *tmax || t < *tmin) )
        return std::nullopt;

    if (t > epsilon)
        // return ray_origin + ray_vector * t;
        return t;

    return std::nullopt;
}