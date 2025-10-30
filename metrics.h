#ifndef METRICS_H
#define METRICS_H
/*
 * Connected components and clustering metrics.
 * Plain C: int for indices; long long for large pair counts.
 */

#include "graph.h"

/* Component analysis result. */
typedef struct {
    int num_components;   /* number of components found */
    int* comp_id;         /* comp_id[i] = component index for node i */
    int* comp_sizes;      /* comp_sizes[c] = size of component c */
    int giant_size;       /* size of largest component */
    long long total_pairs;/* sum over components of C(size,2) */
} Components;

void components_init(Components* comp, int n);
void components_free(Components* comp);

/* Find connected components via BFS over the graph. */
void find_components(Graph* g, Components* comp);

/*
 * Average LOCAL clustering coefficient (Watts-Strogatz definition).
 * For node i with degree k_i and e_i edges among its neighbors:
 *   C_i = 2*e_i / (k_i*(k_i-1))  if k_i >= 2, else C_i = 0
 * Returns: C = (1/N) * sum_i(C_i) over all N nodes.
 */
double compute_clustering_local(Graph* g);

/*
 * Global clustering coefficient (transitivity).
 * C_global = (3 * #triangles) / #wedges
 * where #wedges = sum_i k_i*(k_i-1)/2
 */
double compute_clustering_global(Graph* g);

/*
 * Default clustering function (calls compute_clustering_local by default).
 * This is the function called by sim.c and exported to WASM.
 */
double compute_clustering(Graph* g);

#endif /* METRICS_H */
