#ifndef VEC3_H
#define VEC3_H

#include <math.h>

typedef struct {
    double x, y, z;
} Vec3;

// Constructor
static inline Vec3 vec3(double x, double y, double z) {
    return (Vec3){x, y, z};
}

// Basic operations
static inline Vec3 vec3_add(Vec3 a, Vec3 b) {
    return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline Vec3 vec3_mul(Vec3 a, Vec3 b) {
    return (Vec3){a.x * b.x, a.y * b.y, a.z * b.z};
}

static inline Vec3 vec3_scale(Vec3 v, double t) {
    return (Vec3){v.x * t, v.y * t, v.z * t};
}

static inline Vec3 vec3_negate(Vec3 v) {
    return (Vec3){-v.x, -v.y, -v.z};
}

// Dot and cross products
static inline double vec3_dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline Vec3 vec3_cross(Vec3 a, Vec3 b) {
    return (Vec3){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

// Length operations
static inline double vec3_length_squared(Vec3 v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

static inline double vec3_length(Vec3 v) {
    return sqrt(vec3_length_squared(v));
}

static inline Vec3 vec3_normalize(Vec3 v) {
    double len = vec3_length(v);
    if (len > 0.0) {
        return vec3_scale(v, 1.0 / len);
    }
    return v;
}

// Utility
static inline double vec3_distance(Vec3 a, Vec3 b) {
    return vec3_length(vec3_sub(a, b));
}

// Reflect vector v around normal n
static inline Vec3 vec3_reflect(Vec3 v, Vec3 n) {
    return vec3_sub(v, vec3_scale(n, 2.0 * vec3_dot(v, n)));
}

// Linear interpolation
static inline Vec3 vec3_lerp(Vec3 a, Vec3 b, double t) {
    return vec3_add(vec3_scale(a, 1.0 - t), vec3_scale(b, t));
}

#endif // VEC3_H
