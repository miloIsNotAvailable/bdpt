#ifndef STRUCTS_H
#define STRUCTS_H

#include <cmath>

struct vec3 
{
    float x, y, z;
    vec3(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}

    vec3 operator+(const vec3& other) const {
        return vec3(x + other.x, y + other.y, z + other.z);
    }

    vec3 operator-( const vec3 &b ) const {
        return vec3( this->x - b.x, this->y - b.y, this->z - b.z );
    }

    vec3 operator*( const float &t ) const {
        return vec3( x * t, y * t, z * t );
    }

    static float dot(const vec3& a, const vec3& b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    static vec3 cross(const vec3& a, const vec3& b) {
        return vec3(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        );
    }

    void normalize() {
        float len = std::sqrt(x*x + y*y + z*z);
        if (len > 0.00001f) {
            x /= len; y /= len; z /= len;
        }
    }
};

struct triangle3 {
    vec3 a, b, c;

    triangle3( vec3 a, vec3 b, vec3 c ): a(a), b(b), c(c) {}

    triangle3 operator-( const triangle3 &b ) const {
        return triangle3( this->a - b.a, this->b - b.b, this->c - b.c );
    }

};

inline void perspective(float fov, float aspect, float near, float far, float* mat) {
    float f = 1.0f / std::tan(fov * 0.5f);

    for (int i = 0; i < 16; i++) mat[i] = 0.0f;

    mat[0] = f / aspect;
    mat[5] = f;
    mat[10] = (far + near) / (near - far);
    mat[11] = -1.0f;
    mat[14] = (2.0f * far * near) / (near - far);
}

// Create a lookAt view matrix (column-major order)
inline void lookAt(const vec3& eye, const vec3& center, const vec3& up, float* mat) {
    vec3 f = center - eye;
    f.normalize();

    vec3 upN = up;
    upN.normalize();

    vec3 s = vec3::cross(f, upN);
    s.normalize();

    vec3 u = vec3::cross(s, f);

    for (int i = 0; i < 16; i++) mat[i] = 0.0f;

    mat[0] = s.x;  mat[4] = s.y;  mat[8]  = s.z;
    mat[1] = u.x;  mat[5] = u.y;  mat[9]  = u.z;
    mat[2] = -f.x; mat[6] = -f.y; mat[10] = -f.z;
    mat[15] = 1.0f;

    mat[12] = -vec3::dot(s, eye);
    mat[13] = -vec3::dot(u, eye);
    mat[14] = vec3::dot(f, eye);
}

#endif