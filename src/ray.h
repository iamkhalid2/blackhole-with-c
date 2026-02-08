#ifndef RAY_H
#define RAY_H

#include "vec3.h"

typedef struct {
    Vec3 origin;
    Vec3 direction;  // Should be normalized
} Ray;

// Constructor
static inline Ray ray_create(Vec3 origin, Vec3 direction) {
    return (Ray){origin, vec3_normalize(direction)};
}

// Get point along ray at parameter t
static inline Vec3 ray_at(Ray r, double t) {
    return vec3_add(r.origin, vec3_scale(r.direction, t));
}

#endif // RAY_H
