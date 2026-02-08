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

// Blackbody temperature to RGB color
// Based on Planck spectrum integration over visible wavelengths
// Temperature in Kelvin, returns normalized RGB
static inline Vec3 blackbody_to_rgb(double temperature_K) {
    // Clamp temperature range
    double T = fmax(1000.0, fmin(40000.0, temperature_K));
    
    // Approximate RGB from blackbody (Tanner Helland algorithm, physics-based)
    double r, g, b;
    
    // Normalize temperature to 0-1 range for our working range
    double t = T / 100.0;
    
    // Red channel
    if (t <= 66.0) {
        r = 255.0;
    } else {
        r = 329.698727446 * pow(t - 60.0, -0.1332047592);
    }
    
    // Green channel  
    if (t <= 66.0) {
        g = 99.4708025861 * log(t) - 161.1195681661;
    } else {
        g = 288.1221695283 * pow(t - 60.0, -0.0755148492);
    }
    
    // Blue channel
    if (t >= 66.0) {
        b = 255.0;
    } else if (t <= 19.0) {
        b = 0.0;
    } else {
        b = 138.5177312231 * log(t - 10.0) - 305.0447927307;
    }
    
    // Normalize to 0-1 and clamp
    r = fmax(0.0, fmin(255.0, r)) / 255.0;
    g = fmax(0.0, fmin(255.0, g)) / 255.0;
    b = fmax(0.0, fmin(255.0, b)) / 255.0;
    
    return vec3(r, g, b);
}

// Accretion disk temperature profile (Shakura-Sunyaev thin disk model)
// T(r) ∝ r^(-3/4) × (1 - sqrt(r_in/r))^(1/4)
// Peak temperature typically 10^4 to 10^7 K for stellar/supermassive BH
static inline double disk_temperature(double radius, double r_inner, double T_max) {
    if (radius <= r_inner * 1.01) return T_max;
    
    // Shakura-Sunyaev profile
    double r_ratio = r_inner / radius;
    double profile = pow(radius / r_inner, -0.75) * pow(1.0 - sqrt(r_ratio), 0.25);
    
    return T_max * profile;
}

// Starfield background - Pure Black Space with sharp stars
static inline Vec3 scene_starfield(Vec3 direction) {
    Vec3 dir = vec3_normalize(direction);
    double theta, phi;
    dir_to_spherical(dir, &theta, &phi);
    
    // Higher resolution grid for smaller, sharper stars
    int grid_size = 400;
    int grid_x = (int)((phi + PI) / (2.0 * PI) * grid_size);
    int grid_y = (int)(theta / PI * grid_size);
    
    // Default: PURE BLACK (no nebula haze)
    Vec3 bg = vec3(0.0, 0.0, 0.0);
    
    unsigned int cell_hash = hash2d(grid_x, grid_y);
    double star_prob = hash_to_float(cell_hash);
    
    // Fewer but brighter, sharper stars
    if (star_prob > 0.985) {
        unsigned int brightness_hash = hash(cell_hash + 12345);
        double brightness = hash_to_float(brightness_hash);
        
        // Star color variation (subtle)
        unsigned int color_hash = hash(cell_hash + 67890);
        double color_temp = hash_to_float(color_hash);
        
        Vec3 star_color;
        if (color_temp < 0.2) star_color = vec3(0.8, 0.9, 1.0);      // Blue-ish
        else if (color_temp < 0.8) star_color = vec3(1.0, 1.0, 1.0); // White
        else star_color = vec3(1.0, 0.9, 0.7);                       // Yellow-ish
        
        // Intensity - allow stars to be very bright (bloom source)
        double intensity = 1.0 + brightness * 4.0;
        return vec3_scale(star_color, intensity);
    }
    
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
    
    // Doppler factor - full relativistic calculation
    double doppler = 1.0 / (gamma * (1.0 - v * cos_angle));
    
    // Clamp to avoid extreme values
    doppler = fmax(0.2, fmin(doppler, 4.0));
    
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
// Physics-based: Shakura-Sunyaev temperature profile + blackbody radiation
static inline Vec3 get_disk_color_relativistic(Vec3 point, Vec3 ray_dir, 
                                                double time, double spin_dir,
                                                double gravitational_redshift,
                                                double inner_radius, double outer_radius) {
    double radius = sqrt(point.x * point.x + point.z * point.z);
    
    // Disk boundaries check
    if (radius < inner_radius || radius > outer_radius) {
        return vec3(0, 0, 0);
    }
    
    // PHYSICS-BASED TEMPERATURE: Shakura-Sunyaev profile
    double T_max = 15000.0;  // Peak temp ~15000K for blue-white inner edge
    double T = disk_temperature(radius, inner_radius, T_max);
    
    // Blackbody color from temperature
    Vec3 base_color = blackbody_to_rgb(T);
    
    // Brightness scales as T^4 (Stefan-Boltzmann law)
    double T_ratio = T / T_max;
    double base_brightness = T_ratio * T_ratio * T_ratio * T_ratio;
    base_brightness = fmax(0.05, base_brightness);
    // Log exposure - compresses dynamic range for visibility
    base_brightness = log(1.0 + base_brightness * 10.0) / log(11.0);
    base_brightness = fmax(0.2, base_brightness);
    
    // Scale for HDR rendering
    base_color = vec3_scale(base_color, base_brightness * 5.0);
    
    // Subtle spiral arm structure (less pronounced than before)
    double angle = atan2(point.z, point.x);
    double spiral = sin(angle * 2.0 + radius * 0.3 - time * 1.5 * spin_dir);
    double spiral_factor = 0.9 + 0.1 * spiral * spiral;  // Subtle variation
    
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
