#ifndef GRAPH_H
#define GRAPH_H
/*
 * Undirected simple graph with sorted adjacency lists (CSR-like).
 * Plain C version: int types, no const/static/bitwise.
 *
 * Layout:
 *   nbrs[ offs[u] .. offs[u+1]-1 ] is the sorted neighbor list of u.
 * Fields:
 *   n : number of nodes
 *   m : number of undirected edges
 */

#include "rng.h"

typedef struct {
    int n;       /* nodes */
    int m;       /* edges (undirected) */
    int* offs;   /* length n+1 */
    int* nbrs;   /* length 2*m */
} Graph;

/* Allocate basic arrays for a graph of size n (offs length n+1). */
void graph_init(Graph* g, int n);

/* Free all internal arrays. */
void graph_free(Graph* g);

/* Clear to empty (keeps n and offsets array allocated/zeroed). */
void graph_clear(Graph* g);

/*
 * Build a Watts–Strogatz graph:
 *   - Start with a k-regular ring lattice
 *   - Rewire each edge endpoint with probability beta
 * Requirements: n > 0, k is even, 0 <= beta <= 1
 * RNG used for deterministic behavior.
 */
void ws_build(Graph* g, int n, int k, float beta, RNG* rng);

/* Degree of node u (0 if u out of range). */
int graph_degree(Graph* g, int u);

/*
 * Return pointer to u's neighbor list and write its degree to *deg.
 * Returns 0 and *deg=0 if u is out of range.
 */
int* graph_neighbors(Graph* g, int u, int* deg);

/* Test existence of undirected edge (u, v). Returns 1 if present, 0 otherwise. */
int graph_has_edge(Graph* g, int u, int v);

#endif /* GRAPH_H */
