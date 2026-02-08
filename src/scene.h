#ifndef SCENE_H
#define SCENE_H

#include "vec3.h"
#include "ray.h"
#include <math.h>

#define PI 3.14159265358979323846

// Simple hash function for procedural star generation
static inline unsigned int hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

static inline unsigned int hash2d(int x, int y) {
    return hash(x * 374761393 + y * 668265263);
}

static inline double hash_to_float(unsigned int h) {
    return (double)(h & 0xFFFFFF) / (double)0xFFFFFF;
}

// Convert direction to spherical coordinates
static inline void dir_to_spherical(Vec3 dir, double *theta, double *phi) {
    *theta = acos(fmax(-1.0, fmin(1.0, dir.y)));
    *phi = atan2(dir.z, dir.x);
}

// Starfield background
static inline Vec3 scene_starfield(Vec3 direction) {
    Vec3 dir = vec3_normalize(direction);
    double theta, phi;
    dir_to_spherical(dir, &theta, &phi);
    
    int grid_size = 200;
    int grid_x = (int)((phi + PI) / (2.0 * PI) * grid_size);
    int grid_y = (int)(theta / PI * grid_size);
    
    Vec3 bg = vec3(0.0, 0.0, 0.02);
    
    unsigned int cell_hash = hash2d(grid_x, grid_y);
    double star_prob = hash_to_float(cell_hash);
    
    if (star_prob > 0.97) {
        unsigned int brightness_hash = hash(cell_hash + 12345);
        double brightness = 0.5 + 0.5 * hash_to_float(brightness_hash);
        
        unsigned int color_hash = hash(cell_hash + 67890);
        double color_temp = hash_to_float(color_hash);
        
        Vec3 star_color;
        if (color_temp < 0.3) {
            star_color = vec3(0.8, 0.85, 1.0);
        } else if (color_temp < 0.7) {
            star_color = vec3(1.0, 1.0, 1.0);
        } else if (color_temp < 0.9) {
            star_color = vec3(1.0, 0.95, 0.8);
        } else {
            star_color = vec3(1.0, 0.7, 0.4);
        }
        
        return vec3_scale(star_color, brightness * 2.0);
    }
    
    double milky_way = exp(-pow(dir.y * 3.0, 2.0)) * 0.03;
    bg = vec3_add(bg, vec3(milky_way * 0.8, milky_way * 0.7, milky_way * 1.0));
    
    return bg;
}

// Accretion disk with relativistic effects
typedef struct {
    Vec3 center;
    double inner_radius;
    double outer_radius;
    double spin_direction;  // +1 = counterclockwise from above, -1 = clockwise
} AccretionDisk;

// Keplerian orbital velocity at radius r
// v = sqrt(GM/r), in our units with M=1, v = 1/sqrt(r)
static inline double disk_orbital_velocity(double radius) {
    if (radius < 0.1) return 0.0;
    return 1.0 / sqrt(radius);
}

// Calculate Doppler factor for disk element
// δ = 1 / (γ × (1 - β·n)) where β = v/c, n = direction to observer
// Returns brightness multiplier (δ^3 for intensity)
static inline double calculate_doppler_factor(Vec3 disk_point, Vec3 ray_dir, double spin_dir) {
    double radius = sqrt(disk_point.x * disk_point.x + disk_point.z * disk_point.z);
    if (radius < 0.1) return 1.0;
    
    // Velocity direction: perpendicular to radius, in XZ plane
    // For counterclockwise rotation (spin_dir = +1): v = (-z, 0, x) / r
    Vec3 vel_dir = vec3_normalize(vec3(-disk_point.z * spin_dir, 0.0, disk_point.x * spin_dir));
    
    // Orbital velocity magnitude (as fraction of c)
    double v = disk_orbital_velocity(radius);
    v = fmin(v, 0.8);  // Cap at 0.8c for stability
    
    // Lorentz factor
    double gamma = 1.0 / sqrt(1.0 - v * v);
    
    // Angle between velocity and ray direction (ray going TO observer)
    // We want the ray direction pointing TOWARD camera
    Vec3 to_observer = vec3_negate(ray_dir);
    double cos_angle = vec3_dot(vel_dir, to_observer);
    
    // Doppler factor
    double doppler = 1.0 / (gamma * (1.0 - v * cos_angle));
    
    // Clamp to avoid extreme values
    doppler = fmax(0.1, fmin(doppler, 5.0));
    
    return doppler;
}

// Apply gravitational redshift to color
// Shifts color toward red based on how deep in gravitational well
static inline Vec3 apply_redshift(Vec3 color, double redshift_factor) {
    if (redshift_factor >= 1.0) return color;
    if (redshift_factor <= 0.0) return vec3(0, 0, 0);
    
    // Redshift shifts wavelengths: blue -> green -> yellow -> orange -> red
    // We simulate by shifting color channels
    double shift = 1.0 - redshift_factor;  // 0 = no shift, 1 = max shift
    
    // Move blue toward green/yellow, green toward red
    Vec3 shifted;
    shifted.x = color.x + color.y * shift * 0.5 + color.z * shift * 0.3;
    shifted.y = color.y * (1.0 - shift * 0.3) + color.z * shift * 0.2;
    shifted.z = color.z * (1.0 - shift * 0.8);
    
    // Also dim due to energy loss (intensity ~ 1/λ^4 roughly)
    double dim_factor = pow(redshift_factor, 2.0);
    
    return vec3_scale(shifted, dim_factor);
}

// Get accretion disk color with all relativistic effects
static inline Vec3 get_disk_color_relativistic(Vec3 point, Vec3 ray_dir, 
                                                double time, double spin_dir,
                                                double gravitational_redshift) {
    double radius = sqrt(point.x * point.x + point.z * point.z);
    
    // Base temperature color (black body approximation)
    double inner_r = 3.0;
    double outer_r = 12.0;
    double r_norm = (radius - inner_r) / (outer_r - inner_r);
    r_norm = fmax(0.0, fmin(1.0, r_norm));
    
    // Temperature gradient: hot blue-white inside, cooler orange-red outside
    Vec3 inner_color = vec3(1.5, 1.8, 3.0);
    Vec3 mid_color = vec3(2.5, 1.5, 0.5);
    Vec3 outer_color = vec3(1.5, 0.3, 0.1);
    
    Vec3 base_color;
    if (r_norm < 0.5) {
        base_color = vec3_lerp(inner_color, mid_color, r_norm * 2.0);
    } else {
        base_color = vec3_lerp(mid_color, outer_color, (r_norm - 0.5) * 2.0);
    }
    
    // Base brightness (T^4 scaling approximately)
    double temp_factor = 1.0 - r_norm * 0.8;
    double base_brightness = temp_factor * temp_factor * temp_factor * temp_factor;
    base_brightness = fmax(0.1, base_brightness);
    
    // Spiral arm structure
    double angle = atan2(point.z, point.x);
    double spiral = sin(angle * 3.0 + radius * 0.5 - time * 2.0 * spin_dir);
    double spiral_factor = 0.7 + 0.3 * spiral * spiral;
    
    // Turbulence noise
    unsigned int noise_x = (unsigned int)((point.x + 100.0) * 10.0);
    unsigned int noise_z = (unsigned int)((point.z + 100.0) * 10.0);
    unsigned int noise_hash = hash2d(noise_x, noise_z);
    double noise = 0.8 + 0.4 * hash_to_float(noise_hash);
    
    // Apply base effects
    Vec3 color = vec3_scale(base_color, base_brightness * spiral_factor * noise);
    
    // RELATIVISTIC EFFECTS:
    
    // 1. Doppler beaming (approaching side brighter)
    double doppler = calculate_doppler_factor(point, ray_dir, spin_dir);
    double doppler_brightness = doppler * doppler * doppler;  // I' = I * δ^3
    
    // 2. Doppler color shift (approaching = blueshift, receding = redshift)
    if (doppler > 1.0) {
        // Blueshifted - shift toward blue/white
        double blue_shift = (doppler - 1.0) * 0.3;
        color.z += color.x * blue_shift;
        color.y += color.x * blue_shift * 0.5;
    } else {
        // Redshifted - shift toward red
        double red_shift = (1.0 - doppler) * 0.5;
        color.x += (color.y + color.z) * red_shift;
        color.y *= (1.0 - red_shift * 0.5);
        color.z *= (1.0 - red_shift * 0.8);
    }
    
    // Apply Doppler brightness
    color = vec3_scale(color, doppler_brightness);
    
    // 3. Gravitational redshift
    color = apply_redshift(color, gravitational_redshift);
    
    return color;
}

// Check if ray crosses disk plane
static inline int check_disk_crossing(Vec3 pos1, Vec3 pos2, Vec3 *intersection) {
    if ((pos1.y > 0 && pos2.y <= 0) || (pos1.y < 0 && pos2.y >= 0)) {
        double t = -pos1.y / (pos2.y - pos1.y);
        *intersection = vec3_lerp(pos1, pos2, t);
        return 1;
    }
    return 0;
}

#endif // SCENE_H
