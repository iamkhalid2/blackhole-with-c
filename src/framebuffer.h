#ifndef FRAMEBUFFER_H
#define FRAMEBUFFER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec3.h"

typedef struct {
    int width;
    int height;
    Vec3 *pixels;  // RGB in [0,1] range (HDR, will tonemap on output)
} Framebuffer;

// Create framebuffer
static inline Framebuffer *framebuffer_create(int width, int height) {
    Framebuffer *fb = (Framebuffer *)malloc(sizeof(Framebuffer));
    if (!fb) return NULL;
    
    fb->width = width;
    fb->height = height;
    fb->pixels = (Vec3 *)calloc(width * height, sizeof(Vec3));
    
    if (!fb->pixels) {
        free(fb);
        return NULL;
    }
    return fb;
}

// Free framebuffer
static inline void framebuffer_destroy(Framebuffer *fb) {
    if (fb) {
        free(fb->pixels);
        free(fb);
    }
}

// Set pixel color (HDR values allowed)
static inline void framebuffer_set_pixel(Framebuffer *fb, int x, int y, Vec3 color) {
    if (x >= 0 && x < fb->width && y >= 0 && y < fb->height) {
        fb->pixels[y * fb->width + x] = color;
    }
}

// Get pixel
static inline Vec3 framebuffer_get_pixel(Framebuffer *fb, int x, int y) {
    if (x >= 0 && x < fb->width && y >= 0 && y < fb->height) {
        return fb->pixels[y * fb->width + x];
    }
    return vec3(0, 0, 0);
}

// Clamp value to [0,1]
static inline double clamp01(double v) {
    if (v < 0.0) return 0.0;
    if (v > 1.0) return 1.0;
    return v;
}

// Gamma correction (gamma = 2.2)
static inline double gamma_correct(double v) {
    return pow(clamp01(v), 1.0 / 2.2);
}

// ACES tonemapping (simple approximation)
static inline double tonemap_aces(double x) {
    // Attempt ACES filmic curve
    double a = 2.51;
    double b = 0.03;
    double c = 2.43;
    double d = 0.59;
    double e = 0.14;
    return clamp01((x * (a * x + b)) / (x * (c * x + d) + e));
}

// Write framebuffer to PPM file
static inline int framebuffer_write_ppm(Framebuffer *fb, const char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "Error: Could not open file %s for writing\n", filename);
        return 0;
    }
    
    // PPM header
    fprintf(f, "P3\n%d %d\n255\n", fb->width, fb->height);
    
    // Pixel data (top to bottom)
    for (int y = 0; y < fb->height; y++) {
        for (int x = 0; x < fb->width; x++) {
            Vec3 c = framebuffer_get_pixel(fb, x, y);
            
            // Tonemap and gamma correct
            int r = (int)(255.0 * gamma_correct(tonemap_aces(c.x)));
            int g = (int)(255.0 * gamma_correct(tonemap_aces(c.y)));
            int b = (int)(255.0 * gamma_correct(tonemap_aces(c.z)));
            
            fprintf(f, "%d %d %d\n", r, g, b);
        }
    }
    
    fclose(f);
    return 1;
}

#endif // FRAMEBUFFER_H
