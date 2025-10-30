#ifndef LAYOUT_H
#define LAYOUT_H
/*
 * Layout engine providing:
 *   - Ring layout (nodes on a circle)
 *   - Simple bounded force-directed steps (Fruchterman–Reingold-ish)
 *
 * Plain C:
 *   - int for counts and indices
 *   - int type for layout type (0=RING, 1=FORCE)
 *   - no const/static/bitwise
 */

#include "graph.h"
#include "rng.h"

/* Layout type constants (use ints rather than enum). */
#define LAYOUT_RING  0
#define LAYOUT_FORCE 1

typedef struct {
    int n;
    int type;
    float width, height;

    float* pos_x;
    float* pos_y;
    float* vel_x;
    float* vel_y;

    /* Force parameters */
    float k_rep;          /* repulsion strength */
    float k_spring;       /* spring strength along edges */
    float damping;        /* velocity damping per step */
    float dt;             /* time step */
    float desired_length; /* target edge length */
} Layout;

/* Allocate and initialize layout buffers for n nodes. */
void layout_init(Layout* layout, int n, int type, float width, float height);

/* Free all internal arrays. */
void layout_free(Layout* layout);

/* Update the viewport size and scale parameters accordingly. */
void layout_set_viewport(Layout* layout, float width, float height);

/* Place nodes evenly on a circle centered in the viewport. */
void layout_compute_ring(Layout* layout);

/* Initialize a compact randomized placement for force layout. */
void layout_init_force(Layout* layout, RNG* rng);

/* Execute one simulation step of force-directed layout (bounded). */
void layout_step_force(Layout* layout, Graph* g);

#endif /* LAYOUT_H */
