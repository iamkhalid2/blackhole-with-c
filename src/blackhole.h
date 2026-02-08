#ifndef BLACKHOLE_H
#define BLACKHOLE_H

#include "vec3.h"
#include "ray.h"
#include <math.h>

#define PI 3.14159265358979323846

// Kerr black hole with spin
typedef struct {
    Vec3 position;
    double mass;
    double spin;                    // Dimensionless spin parameter a/M, range [0, 1)
    double schwarzschild_radius;    // Rs = 2M
    double event_horizon_radius;    // r+ = M + sqrt(M^2 - a^2)
    Vec3 spin_axis;                 // Direction of angular momentum (usually +Y)
} BlackHole;

// Create Kerr black hole
static inline BlackHole blackhole_create(Vec3 position, double mass, double spin) {
    BlackHole bh;
    bh.position = position;
    bh.mass = mass;
    bh.spin = fmin(0.998, fmax(0.0, spin));  // Clamp to avoid singularities
    bh.schwarzschild_radius = 2.0 * mass;
    
    // Kerr event horizon: r+ = M + sqrt(M^2 - a^2*M^2) = M(1 + sqrt(1 - a^2))
    double a = bh.spin;
    bh.event_horizon_radius = mass * (1.0 + sqrt(1.0 - a * a));
    
    bh.spin_axis = vec3(0, 1, 0);  // Spin around Y axis
    return bh;
}

// Calculate Innermost Stable Circular Orbit (ISCO) radius for Kerr black hole
// This is where the accretion disk inner edge physically must be
// Formula from Bardeen, Press, Teukolsky (1972)
static inline double blackhole_isco_radius(BlackHole *bh) {
    double a = bh->spin;
    double M = bh->mass;
    
    // Z1 = 1 + (1-a²)^(1/3) × [(1+a)^(1/3) + (1-a)^(1/3)]
    double one_minus_a2 = 1.0 - a * a;
    double cbrt_1ma2 = cbrt(one_minus_a2);
    double cbrt_1pa = cbrt(1.0 + a);
    double cbrt_1ma = cbrt(1.0 - a);
    double Z1 = 1.0 + cbrt_1ma2 * (cbrt_1pa + cbrt_1ma);
    
    // Z2 = sqrt(3a² + Z1²)
    double Z2 = sqrt(3.0 * a * a + Z1 * Z1);
    
    // r_ISCO = M × (3 + Z2 - sqrt((3-Z1)(3+Z1+2Z2)))
    // This is for prograde orbits (co-rotating with BH)
    double r_isco = M * (3.0 + Z2 - sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    
    return r_isco;
}

// Frame dragging angular velocity (Lense-Thirring effect)
// At distance r, space itself rotates with angular velocity omega
static inline double blackhole_frame_drag_omega(BlackHole *bh, double r, double theta) {
    if (r <= bh->event_horizon_radius) return 0.0;
    
    double a = bh->spin * bh->mass;
    double r2 = r * r;
    double a2 = a * a;
    double sin_theta = sin(theta);
    double sin2 = sin_theta * sin_theta;
    
    // Simplified frame dragging: omega = 2Mar / (r^2 + a^2)^2
    // This is the leading order term
    double denom = (r2 + a2) * (r2 + a2);
    if (denom < 1e-10) return 0.0;
    
    return 2.0 * bh->mass * a * r / denom;
}

// Calculate gravitational redshift factor at radius r
// z = 1/sqrt(1 - Rs/r) for Schwarzschild
// For Kerr, it's more complex but we use approximation
static inline double blackhole_redshift_factor(BlackHole *bh, double r) {
    if (r <= bh->event_horizon_radius) return 0.0;
    
    double rs = bh->schwarzschild_radius;
    double factor = 1.0 - rs / r;
    if (factor <= 0.0) return 0.0;
    
    return sqrt(factor);
}

// Get polar angle theta from position (angle from spin axis)
static inline double get_theta(BlackHole *bh, Vec3 pos) {
    Vec3 rel = vec3_sub(pos, bh->position);
    double r = vec3_length(rel);
    if (r < 1e-10) return PI / 2.0;
    
    // Theta is angle from spin axis (Y)
    double cos_theta = vec3_dot(vec3_normalize(rel), bh->spin_axis);
    return acos(fmax(-1.0, fmin(1.0, cos_theta)));
}

// Bend ray direction with Kerr corrections
// Adds frame dragging twist in addition to gravitational bending
static inline Vec3 blackhole_bend_ray_kerr(BlackHole *bh, Vec3 ray_pos, Vec3 ray_dir, double step_size) {
    Vec3 to_bh = vec3_sub(bh->position, ray_pos);
    double dist = vec3_length(to_bh);
    
    if (dist < bh->event_horizon_radius * 0.5) {
        return ray_dir;  // Too close, avoid singularity
    }
    
    // Standard gravitational bending (enhanced for Kerr)
    Vec3 gravity_dir = vec3_normalize(to_bh);
    
    // Bending strength - stronger for Kerr near horizon
    double r_ratio = bh->schwarzschild_radius / dist;
    double bending_strength = 1.5 * r_ratio * r_ratio;
    
    // Kerr enhancement: stronger bending in equatorial plane
    double theta = get_theta(bh, ray_pos);
    double equatorial_boost = 1.0 + 0.5 * bh->spin * sin(theta) * sin(theta);
    bending_strength *= equatorial_boost;
    
    // Perpendicular component for bending
    double parallel = vec3_dot(gravity_dir, ray_dir);
    Vec3 perp = vec3_sub(gravity_dir, vec3_scale(ray_dir, parallel));
    
    Vec3 new_dir = vec3_add(ray_dir, vec3_scale(perp, bending_strength * step_size));
    
    // Frame dragging: twist ray direction around spin axis
    double omega = blackhole_frame_drag_omega(bh, dist, theta);
    if (omega > 1e-10) {
        // Rotate ray direction around spin axis by omega * step_size
        double angle = omega * step_size;
        double cos_a = cos(angle);
        double sin_a = sin(angle);
        
        Vec3 axis = bh->spin_axis;
        
        // Rodrigues rotation formula
        Vec3 rotated = vec3_add(
            vec3_add(
                vec3_scale(new_dir, cos_a),
                vec3_scale(vec3_cross(axis, new_dir), sin_a)
            ),
            vec3_scale(axis, vec3_dot(axis, new_dir) * (1.0 - cos_a))
        );
        new_dir = rotated;
    }
    
    return vec3_normalize(new_dir);
}

// Check if position is inside event horizon
// For Kerr, horizon shape depends on theta
static inline int blackhole_is_inside_horizon(BlackHole *bh, Vec3 pos) {
    Vec3 rel = vec3_sub(pos, bh->position);
    double r = vec3_length(rel);
    
    // Kerr horizon is oblate - larger at equator
    double theta = get_theta(bh, pos);
    double a = bh->spin * bh->mass;
    
    // r+ = M + sqrt(M^2 - a^2 cos^2(theta)) approximately
    // For visual purposes, use simpler oblate shape
    double cos_theta = cos(theta);
    double horizon_r = bh->event_horizon_radius * (1.0 - 0.3 * bh->spin * cos_theta * cos_theta);
    
    return r <= horizon_r;
}

// Result of tracing a ray near black hole
typedef struct {
    int swallowed;
    Vec3 final_direction;
    Vec3 final_position;
    double total_distance;
    double min_radius;        // Closest approach for redshift
} BlackHoleTraceResult;

// Trace ray through curved spacetime (Kerr)
static inline BlackHoleTraceResult blackhole_trace_ray(BlackHole *bh, Ray ray, double max_distance) {
    BlackHoleTraceResult result = {0, ray.direction, ray.origin, 0.0, 1e10};
    
    Vec3 pos = ray.origin;
    Vec3 dir = ray.direction;
    double total_dist = 0.0;
    double min_r = 1e10;
    
    double base_step = 0.05;
    
    while (total_dist < max_distance) {
        double dist_to_bh = vec3_distance(pos, bh->position);
        min_r = fmin(min_r, dist_to_bh);
        
        // Check if swallowed
        if (blackhole_is_inside_horizon(bh, pos)) {
            result.swallowed = 1;
            result.total_distance = total_dist;
            result.min_radius = min_r;
            return result;
        }
        
        // Adaptive step size
        double step = base_step * fmax(0.01, dist_to_bh / (bh->schwarzschild_radius * 5.0));
        step = fmin(step, 0.5);
        
        // Bend ray with Kerr physics
        dir = blackhole_bend_ray_kerr(bh, pos, dir, step);
        
        // Move along ray
        pos = vec3_add(pos, vec3_scale(dir, step));
        total_dist += step;
        
        // Check if escaped
        if (dist_to_bh > bh->schwarzschild_radius * 25.0) {
            Vec3 to_bh = vec3_sub(bh->position, pos);
            if (vec3_dot(dir, to_bh) < 0) {
                break;
            }
        }
    }
    
    result.final_direction = dir;
    result.final_position = pos;
    result.total_distance = total_dist;
    result.min_radius = min_r;
    return result;
}

#endif // BLACKHOLE_H
