#ifndef WS_H
#define WS_H
/*
 * Watts–Strogatz Small-World Graph (Plain C / WASM-friendly)
 *
 * Design goals for this header:
 * - Plain C types only: use `int` and pointers (no const/unsigned/bool).
 * - Minimal includes here; implementation includes stdlib/string.
 * - Keep a small, clear public API for building a ring lattice, rewiring,
 *   exporting edges, and basic accessors.
 *
 * Graph Model:
 *   Undirected, simple graph on n nodes labeled [0, n-1].
 *   Adjacency is stored as a single neighbors[] array with an offsets[]
 *   index into neighbors for each node:
 *
 *     neighbors[ offsets[u] ... offsets[u+1]-1 ]   => sorted neighbors of u
 *
 *   The array for each node is kept sorted so that we can use binary search
 *   to test edge existence efficiently and keep the representation simple.
 */

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------
 * Public data structures
 * ----------------------------- */

/* Graph storage for an undirected simple graph. */
typedef struct {
    int n;           /* number of nodes (>= 0) */
    int k;           /* even degree used for ring build (each node has degree k) */
    int *offsets;    /* length n+1; offsets[u] is start index into neighbors */
    int *neighbors;  /* length 2E (E = number of undirected edges) */
    int capacity;    /* allocated length of neighbors[] (in ints) */
} WSGraph;

/* Very simple linear-congruential PRNG state. */
typedef struct {
    int state;       /* non-zero seed */
} WSPrng;

/* -----------------------------
 * PRNG interface (plain arithmetic; no bitwise)
 * ----------------------------- */

/* Seed the PRNG (if seed==0, it is promoted to 1 to avoid the zero cycle). */
void ws_prng_seed(WSPrng *rng, int seed);

/* Return a pseudo-random float in [0, 1). */
float ws_prng_next(WSPrng *rng);

/* Return a pseudo-random integer in [0, max). If max==0, returns 0. */
int ws_prng_nextInt(WSPrng *rng, int max);

/* -----------------------------
 * Graph construction and mutation
 * ----------------------------- */

/*
 * Build a k-regular ring lattice on n nodes.
 * Each node u is connected to k/2 neighbors on each side around the ring.
 * Requirements: n > 0, k is even, and k < n.
 * On success, returns 0 and fills *out (caller must call ws_free).
 * On error, returns <0 and leaves *out in a zeroed state.
 */
int ws_build_ring(int n, int k, WSGraph *out);

/*
 * Perform a single rewiring pass:
 * Rewires approximately (frac * beta * E) edges chosen uniformly.
 *   - Avoids self-loops and duplicate edges.
 *   - Preserves undirected symmetry and sorted neighbor lists.
 * Returns: number of edges rewired (>=0), or <0 on error.
 */
int ws_rewire_pass(WSGraph *g, WSPrng *rng, float beta, float frac);

/*
 * Export the unique undirected edge list as pairs (u, v) with u < v.
 * The out_pairs buffer must have room for at least E pairs (2*E ints).
 * Returns: number of edges E on success, or:
 *   -1 if g or out_pairs is NULL
 *   -2 if max_pairs < E
 */
int ws_export_edgelist(WSGraph *g, int *out_pairs, int max_pairs);

/* -----------------------------
 * Basic accessors
 * ----------------------------- */

/* Number of nodes and edges (edges are undirected; E = |neighbors|/2). */
int ws_num_nodes(WSGraph *g);
int ws_num_edges(WSGraph *g);

/* Zero-copy access to CSR-like arrays. (Non-const to keep API plain-C.) */
int* ws_offsets(WSGraph *g);
int* ws_neighbors(WSGraph *g);

/* Free internal arrays (does not free the WSGraph struct itself). */
void ws_free(WSGraph *g);

/* -----------------------------
 * Metrics computation
 * ----------------------------- */

/*
 * Compute average clustering coefficient (Watts-Strogatz formula).
 *
 * For each node i with degree k_i >= 2:
 *   C_i = 2 * e_i / (k_i * (k_i - 1))
 * where e_i is the number of edges among i's neighbors (triangles at i).
 *
 * Returns: C = (1/N) * sum(C_i) for nodes with degree >= 2.
 * Returns 0.0 if graph is NULL or has no nodes with degree >= 2.
 *
 * Algorithm: Two-pointer set intersection on sorted neighbor arrays.
 * Time complexity: O(sum of d_i^2) where d_i is degree of node i.
 * Space complexity: O(1) auxiliary space.
 */
double ws_avg_clustering(WSGraph *g);

#ifdef __cplusplus
}
#endif
#endif /* WS_H */
