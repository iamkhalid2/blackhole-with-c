#ifndef CAMERA_H
#define CAMERA_H

#include "vec3.h"
#include "ray.h"
#include <math.h>

typedef struct {
    Vec3 position;
    Vec3 forward;
    Vec3 right;
    Vec3 up;
    double fov_radians;
    double aspect_ratio;
} Camera;

// Create camera looking at target
static inline Camera camera_create(Vec3 position, Vec3 look_at, Vec3 world_up, double fov_degrees, double aspect_ratio) {
    Camera cam;
    cam.position = position;
    cam.fov_radians = fov_degrees * 3.14159265358979323846 / 180.0;
    cam.aspect_ratio = aspect_ratio;
    
    // Compute basis vectors
    cam.forward = vec3_normalize(vec3_sub(look_at, position));
    cam.right = vec3_normalize(vec3_cross(cam.forward, world_up));
    cam.up = vec3_cross(cam.right, cam.forward);
    
    return cam;
}

// Get ray for normalized screen coordinates (u, v in [0, 1])
static inline Ray camera_get_ray(Camera *cam, double u, double v) {
    // Convert to [-1, 1] range, centered
    double x = (2.0 * u - 1.0) * cam->aspect_ratio;
    double y = (2.0 * v - 1.0);
    
    // Scale by FOV
    double scale = tan(cam->fov_radians / 2.0);
    x *= scale;
    y *= scale;
    
    // Compute ray direction in world space
    Vec3 dir = vec3_add(
        vec3_add(
            cam->forward,
            vec3_scale(cam->right, x)
        ),
        vec3_scale(cam->up, y)
    );
    
    return ray_create(cam->position, dir);
}

#endif // CAMERA_H
