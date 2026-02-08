/*
 * Black Hole Simulation in Pure C
 * 
 * Kerr (rotating) black hole with:
 * - Frame dragging (Lense-Thirring effect)
 * - Relativistic Doppler beaming
 * - Gravitational redshift
 * - Asymmetric (D-shaped) shadow
 * 
 * Output: PPM image files
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "vec3.h"
#include "ray.h"
#include "camera.h"
#include "framebuffer.h"
#include "blackhole.h"
#include "scene.h"

// Configuration
#define TEST_MODE 0         // 1 = low-res fast test, 0 = full quality

#if TEST_MODE
#define WIDTH 200
#define HEIGHT 150
#else
#define WIDTH 1600
#define HEIGHT 900
#endif

#define MAX_RAY_DISTANCE 150.0

// Animation settings
#define ANIMATE 1           // Set to 1 for animation, 0 for single frame
#define NUM_FRAMES 30       // Render 30 frames as requested
#define ORBIT_RADIUS 120.0  // Very far back - full system visible
#define CAMERA_HEIGHT 5.0   // Minimal tilt - near edge-on like Interstellar

// Black hole settings
#define BH_SPIN 0.95        // High spin for visible effects (0 to 0.998)
#define DISK_SPIN_DIR 1.0   // +1 = counterclockwise, -1 = clockwise

// Accretion disk settings
// Note: DISK_INNER_RADIUS will be overridden by physics-based ISCO
#define DISK_INNER_RADIUS 3.0  // Fallback - actual inner edge uses ISCO
#define DISK_OUTER_RADIUS 50.0  // Wide disk spanning frame

void print_progress(int current, int total) {
    int bar_width = 50;
    float progress = (float)current / total;
    int filled = (int)(bar_width * progress);
    
    printf("\r[");
    for (int i = 0; i < bar_width; i++) {
        if (i < filled) printf("=");
        else if (i == filled) printf(">");
        else printf(" ");
    }
    printf("] %3d%%", (int)(progress * 100));
    fflush(stdout);
}

// Trace result with disk info
typedef struct {
    int hit_type;           // 0=escaped, 1=swallowed, 2=hit disk
    Vec3 final_direction;
    Vec3 disk_color;
    double min_radius;
} FullTraceResult;

FullTraceResult trace_ray_full(BlackHole *bh, Ray ray, double time) {
    FullTraceResult result = {0, ray.direction, vec3(0,0,0), 1e10};
    
    Vec3 pos = ray.origin;
    Vec3 dir = ray.direction;
    double total_dist = 0.0;
    Vec3 accumulated_disk = vec3(0, 0, 0);
    double disk_opacity = 0.0;
    double min_r = 1e10;
    
    double base_step = 0.05;
    
    while (total_dist < MAX_RAY_DISTANCE && disk_opacity < 0.98) {
        double dist_to_bh = vec3_distance(pos, bh->position);
        min_r = fmin(min_r, dist_to_bh);
        
        // Check event horizon
        if (blackhole_is_inside_horizon(bh, pos)) {
            result.hit_type = 1;
            result.disk_color = accumulated_disk;
            result.min_radius = min_r;
            return result;
        }
        
        // Adaptive step
        double step = base_step * fmax(0.01, dist_to_bh / (bh->schwarzschild_radius * 5.0));
        step = fmin(step, 0.5);
        
        // Kerr ray bending with frame dragging
        dir = blackhole_bend_ray_kerr(bh, pos, dir, step);
        
        Vec3 new_pos = vec3_add(pos, vec3_scale(dir, step));
        
        // Check disk crossing
        Vec3 disk_point;
        if (check_disk_crossing(pos, new_pos, &disk_point)) {
            double disk_radius = sqrt(disk_point.x * disk_point.x + disk_point.z * disk_point.z);
            
            if (disk_radius >= DISK_INNER_RADIUS && disk_radius <= DISK_OUTER_RADIUS) {
                // Calculate gravitational redshift at this point
                double disk_dist = vec3_distance(disk_point, bh->position);
                double grav_redshift = blackhole_redshift_factor(bh, disk_dist);
                
                // Get disk color with full relativistic effects
                Vec3 disk_emission = get_disk_color_relativistic(
                    disk_point, dir, time, DISK_SPIN_DIR, grav_redshift,
                    DISK_INNER_RADIUS, DISK_OUTER_RADIUS
                );
                
                // Composite
                double alpha = 0.7;
                Vec3 attenuated = vec3_scale(disk_emission, (1.0 - disk_opacity) * alpha);
                accumulated_disk = vec3_add(accumulated_disk, attenuated);
                disk_opacity += (1.0 - disk_opacity) * alpha;
            }
        }
        
        pos = new_pos;
        total_dist += step;
        
        // Escape check
        if (dist_to_bh > bh->schwarzschild_radius * 25.0) {
            Vec3 to_bh = vec3_sub(bh->position, pos);
            if (vec3_dot(dir, to_bh) < 0) {
                break;
            }
        }
    }
    
    result.final_direction = dir;
    result.disk_color = accumulated_disk;
    result.min_radius = min_r;
    
    if (disk_opacity > 0.1) {
        result.hit_type = 2;
    }
    
    return result;
}

void render_frame(int frame_num, int total_frames, double time) {
    char filename[256];
    if (total_frames > 1) {
        snprintf(filename, sizeof(filename), "output/frame_%04d.ppm", frame_num);
    } else {
        snprintf(filename, sizeof(filename), "output/blackhole.ppm");
    }
    
    Framebuffer *fb = framebuffer_create(WIDTH, HEIGHT);
    if (!fb) {
        fprintf(stderr, "Error: Could not allocate framebuffer\n");
        return;
    }
    
    // Kerr black hole with spin
    BlackHole bh = blackhole_create(vec3(0, 0, 0), 1.0, BH_SPIN);
    
    // Camera position
    // Camera position - orbit slowly (100 frames for full 360)
    // We decouple orbit speed from total_frames so 15 frames doesn't mean full super-fast spin
    double angle = (2.0 * 3.14159265358979323846 * frame_num) / 100.0;
    Vec3 cam_pos = vec3(
        ORBIT_RADIUS * sin(angle),
        CAMERA_HEIGHT,
        ORBIT_RADIUS * cos(angle)
    );
    Vec3 look_at = vec3(0, 0, 0);
    Vec3 world_up = vec3(0, 1, 0);
    double fov = 50.0;
    double aspect = (double)WIDTH / (double)HEIGHT;
    
    Camera cam = camera_create(cam_pos, look_at, world_up, fov, aspect);
    
    // OpenMP disabled (needs admin to install) - running single threaded
    for (int y = 0; y < HEIGHT; y++) {
        if (total_frames == 1) {
            print_progress(y, HEIGHT);
        }
        
        for (int x = 0; x < WIDTH; x++) {
            double u = (double)x / (double)(WIDTH - 1);
            double v = 1.0 - (double)y / (double)(HEIGHT - 1);
            
            Ray ray = camera_get_ray(&cam, u, v);
            FullTraceResult trace = trace_ray_full(&bh, ray, time);
            
            Vec3 color;
            
            if (trace.hit_type == 1) {
                // Swallowed by black hole
                color = trace.disk_color;
            } else {
                // Starfield background
                Vec3 bg = scene_starfield(trace.final_direction);
                
                // Blend disk emission
                double disk_intensity = vec3_length(trace.disk_color);
                double bg_factor = 1.0 / (1.0 + disk_intensity * 0.5);
                
                color = vec3_add(vec3_scale(bg, bg_factor), trace.disk_color);
            }
            
            framebuffer_set_pixel(fb, x, y, color);
        }
    }
    
    if (total_frames == 1) {
        print_progress(HEIGHT, HEIGHT);
        printf("\n");
    }
    
    framebuffer_write_ppm(fb, filename);
    framebuffer_destroy(fb);
}

int main(void) {
    printf("Kerr Black Hole Simulation\n");
    printf("Resolution: %d x %d\n", WIDTH, HEIGHT);
    printf("Spin parameter: %.3f\n", BH_SPIN);
    printf("===========================\n\n");
    
    clock_t start_time = clock();
    
    int num_frames = ANIMATE ? NUM_FRAMES : 1;
    
    if (ANIMATE) {
        printf("Rendering %d frames...\n\n", num_frames);
    } else {
        printf("Rendering single frame...\n");
    }
    
    for (int frame = 0; frame < num_frames; frame++) {
        double time = (double)frame / 30.0;
        
        if (ANIMATE) {
            printf("Frame %d/%d ", frame + 1, num_frames);
            fflush(stdout);
        }
        
        render_frame(frame, num_frames, time);
        
        if (ANIMATE) {
            clock_t now = clock();
            double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;
            double per_frame = elapsed / (frame + 1);
            double remaining = per_frame * (num_frames - frame - 1);
            printf("(%.1fs remaining)\n", remaining);
        }
    }
    
    clock_t end_time = clock();
    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    
    printf("\n===========================\n");
    printf("Total render time: %.2f seconds\n", total_time);
    
    if (ANIMATE) {
        printf("\nTo create video:\n");
        printf("  ffmpeg -framerate 30 -i output/frame_%%04d.ppm -c:v libx264 -pix_fmt yuv420p output/blackhole.mp4\n");
    } else {
        printf("Output: output/blackhole.ppm\n");
    }
    
    return 0;
}
