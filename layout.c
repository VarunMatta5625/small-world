#include "layout.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void layout_init(Layout* layout, int n, int type, float width, float height) {
    if (!layout) return;
    layout->n = n;
    layout->type = type;
    layout->width = width;
    layout->height = height;

    layout->pos_x = (float*)calloc((size_t)n, sizeof(float));
    layout->pos_y = (float*)calloc((size_t)n, sizeof(float));
    layout->vel_x = (float*)calloc((size_t)n, sizeof(float));
    layout->vel_y = (float*)calloc((size_t)n, sizeof(float));

    /* Set force parameters for stable, compact layouts. */
    if (n <= 0) n = 1;
    float area = width * height;
    float min_dim = (width < height) ? width : height;

    layout->k_rep = 5.0f * sqrtf(area / (float)n);
    layout->k_spring = 0.8f;
    layout->damping = 0.97f;
    layout->dt = 0.15f;
    layout->desired_length = min_dim * 0.05f;
}

void layout_free(Layout* layout) {
    if (!layout) return;
    if (layout->pos_x) free(layout->pos_x);
    if (layout->pos_y) free(layout->pos_y);
    if (layout->vel_x) free(layout->vel_x);
    if (layout->vel_y) free(layout->vel_y);
    layout->pos_x = 0;
    layout->pos_y = 0;
    layout->vel_x = 0;
    layout->vel_y = 0;
}

void layout_set_viewport(Layout* layout, float width, float height) {
    if (!layout) return;
    layout->width = width;
    layout->height = height;

    if (layout->n <= 0) return;
    float area = width * height;
    float min_dim = (width < height) ? width : height;

    layout->k_rep = 5.0f * sqrtf(area / (float)layout->n);
    layout->desired_length = min_dim * 0.05f;
}

void layout_compute_ring(Layout* layout) {
    if (!layout || layout->n <= 0) return;
    float cx = layout->width / 2.0f;
    float cy = layout->height / 2.0f;
    float min_dim = (layout->width < layout->height) ? layout->width : layout->height;
    float radius = 0.40f * min_dim;

    for (int i = 0; i < layout->n; i++) {
        float angle = 2.0f * (float)M_PI * (float)i / (float)layout->n;
        layout->pos_x[i] = cx + radius * cosf(angle);
        layout->pos_y[i] = cy + radius * sinf(angle);
        layout->vel_x[i] = 0.0f;
        layout->vel_y[i] = 0.0f;
    }
}

void layout_init_force(Layout* layout, RNG* rng) {
    if (!layout) return;
    float cx = layout->width / 2.0f;
    float cy = layout->height / 2.0f;
    float min_dim = (layout->width < layout->height) ? layout->width : layout->height;
    float radius = 0.25f * min_dim;

    for (int i = 0; i < layout->n; i++) {
        float angle = 2.0f * (float)M_PI * (float)i / (float)layout->n;
        float px = (rng_float(rng) - 0.5f) * 5.0f;
        float py = (rng_float(rng) - 0.5f) * 5.0f;
        layout->pos_x[i] = cx + radius * cosf(angle) + px;
        layout->pos_y[i] = cy + radius * sinf(angle) + py;
        layout->vel_x[i] = 0.0f;
        layout->vel_y[i] = 0.0f;
    }
}

void layout_step_force(Layout* layout, Graph* g) {
    if (!layout || !g || layout->n <= 0) return;

    float* fx = (float*)calloc((size_t)layout->n, sizeof(float));
    float* fy = (float*)calloc((size_t)layout->n, sizeof(float));

    float cx = layout->width / 2.0f;
    float cy = layout->height / 2.0f;

    /* Short-range repulsion (to avoid overlap) */
    float cutoff_sq = 100.0f * 100.0f;

    for (int i = 0; i < layout->n; i++) {
        for (int j = i + 1; j < layout->n; j++) {
            float dx = layout->pos_x[i] - layout->pos_x[j];
            float dy = layout->pos_y[i] - layout->pos_y[j];
            float dist_sq = dx * dx + dy * dy;

            if (dist_sq < cutoff_sq && dist_sq > 1.0f) {
                float dist = sqrtf(dist_sq);
                float force = layout->k_rep * layout->k_rep / dist;
                float fdx = force * dx / dist;
                float fdy = force * dy / dist;
                fx[i] += fdx; fy[i] += fdy;
                fx[j] -= fdx; fy[j] -= fdy;
            }
        }
    }

    /* Spring attraction along edges */
    for (int i = 0; i < g->n; i++) {
        int deg = 0;
        int* nbrs = graph_neighbors(g, i, &deg);
        for (int k = 0; k < deg; k++) {
            int j = nbrs[k];
            if (i >= j) continue;

            float dx = layout->pos_x[j] - layout->pos_x[i];
            float dy = layout->pos_y[j] - layout->pos_y[i];
            float dist = sqrtf(dx * dx + dy * dy);
            if (dist < 0.1f) dist = 0.1f;

            float force = layout->k_spring * (dist - layout->desired_length);
            float fdx = force * dx / dist;
            float fdy = force * dy / dist;
            fx[i] += fdx; fy[i] += fdy;
            fx[j] -= fdx; fy[j] -= fdy;
        }
    }

    /* Strong centering force toward viewport center */
    float k_center = 0.15f;
    for (int i = 0; i < layout->n; i++) {
        float dx = cx - layout->pos_x[i];
        float dy = cy - layout->pos_y[i];
        float d = sqrtf(dx * dx + dy * dy);
        float cf = k_center * d;
        if (d > 0.1f) {
            fx[i] += cf * dx / d;
            fy[i] += cf * dy / d;
        }
    }

    /* Integrate velocities/positions with soft bounding. */
    float margin = 100.0f;
    float max_vel = 5.0f;
    float bounce_damping = 0.2f;

    for (int i = 0; i < layout->n; i++) {
        layout->vel_x[i] = (layout->vel_x[i] + fx[i] * layout->dt) * layout->damping;
        layout->vel_y[i] = (layout->vel_y[i] + fy[i] * layout->dt) * layout->damping;

        float vm = sqrtf(layout->vel_x[i] * layout->vel_x[i] +
                         layout->vel_y[i] * layout->vel_y[i]);
        if (vm > max_vel) {
            layout->vel_x[i] = layout->vel_x[i] * max_vel / vm;
            layout->vel_y[i] = layout->vel_y[i] * max_vel / vm;
        }

        layout->pos_x[i] += layout->vel_x[i];
        layout->pos_y[i] += layout->vel_y[i];

        if (layout->pos_x[i] < margin) {
            layout->pos_x[i] = margin;
            layout->vel_x[i] = -layout->vel_x[i] * bounce_damping;
        }
        if (layout->pos_x[i] > layout->width - margin) {
            layout->pos_x[i] = layout->width - margin;
            layout->vel_x[i] = -layout->vel_x[i] * bounce_damping;
        }
        if (layout->pos_y[i] < margin) {
            layout->pos_y[i] = margin;
            layout->vel_y[i] = -layout->vel_y[i] * bounce_damping;
        }
        if (layout->pos_y[i] > layout->height - margin) {
            layout->pos_y[i] = layout->height - margin;
            layout->vel_y[i] = -layout->vel_y[i] * bounce_damping;
        }
    }

    free(fx);
    free(fy);
}
