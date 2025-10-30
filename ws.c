/*
 * Watts–Strogatz Small-World Graph (Plain C Implementation)
 *
 * This file intentionally avoids:
 *   - unsigned types (we use plain `int`)
 *   - const and static keywords
 *   - bitwise operations
 *
 * It uses simple, readable algorithms:
 *   - Linear-congruential generator (LCG) for PRNG
 *   - Insertion sort for tiny per-node neighbor lists
 *   - Binary search for edge existence in sorted neighbor lists
 *   - Simple array growth by doubling for neighbors[]
 *
 * Representation invariant:
 *   - neighbors[ offsets[u] .. offsets[u+1]-1 ] is sorted for each u
 *   - Graph is undirected and simple (no duplicates, no self-loops)
 */

#include "ws.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================
 * 1) Utility helpers (small and readable)
 * ============================================================ */

/* Return the number of edges (undirected) based on total neighbor slots. */
int ws_num_edges(WSGraph *g) {
    if (g == 0 || g->offsets == 0) return 0;
    return g->offsets[g->n] / 2;
}

int ws_num_nodes(WSGraph *g) {
    if (g == 0) return 0;
    return g->n;
}

int* ws_offsets(WSGraph *g)   { return g ? g->offsets   : 0; }
int* ws_neighbors(WSGraph *g) { return g ? g->neighbors : 0; }

/* Safe reallocation of neighbors[] to at least need_total ints. */
void ws_ensure_neighbor_capacity(WSGraph *g, int need_total) {
    if (g->capacity >= need_total) return;

    /* Double strategy for simplicity and amortized O(1) growth. */
    if (g->capacity <= 0) g->capacity = 16;
    while (g->capacity < need_total) g->capacity = g->capacity * 2;

    g->neighbors = (int*)realloc(g->neighbors, (size_t)g->capacity * sizeof(int));
}

/* Binary search: does sorted array a[0..len-1] contain x ? (returns 1/0) */
int ws_bs_has(int *a, int len, int x) {
    int lo = 0;
    int hi = len;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        int v = a[mid];
        if (v == x) return 1;
        if (v < x)  lo = mid + 1; else hi = mid;
    }
    return 0;
}

/* Where should x be inserted in sorted array a[0..len-1]? */
int ws_bs_insert_pos(int *a, int len, int x) {
    int lo = 0;
    int hi = len;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (a[mid] < x) lo = mid + 1; else hi = mid;
    }
    return lo;
}

/* Insert neighbor v into u's sorted neighbor list and shift tails. */
void ws_insert_neighbor(WSGraph *g, int u, int v) {
    int start = g->offsets[u];
    int end   = g->offsets[u+1];
    int deg   = end - start;

    int ins_rel = ws_bs_insert_pos(g->neighbors + start, deg, v);
    int ins_abs = start + ins_rel;

    /* Total used slots across all nodes: offsets[n] */
    int total = g->offsets[g->n];
    ws_ensure_neighbor_capacity(g, total + 1);

    /* Move the suffix right by one to make a gap. */
    memmove(g->neighbors + ins_abs + 1,
            g->neighbors + ins_abs,
            (size_t)(total - ins_abs) * sizeof(int));

    /* Write new neighbor at the gap. */
    g->neighbors[ins_abs] = v;

    /* Shift all subsequent offsets by +1. */
    {
        int i;
        for (i = u + 1; i <= g->n; ++i) g->offsets[i] = g->offsets[i] + 1;
    }
}

/* Remove neighbor v from u's list (assumes it exists), then shift tails. */
void ws_erase_neighbor(WSGraph *g, int u, int v) {
    int start = g->offsets[u];
    int end   = g->offsets[u+1];
    int total = g->offsets[g->n];

    /* Find the position linearly (deg is usually small). */
    {
        int i;
        for (i = start; i < end; ++i) {
            if (g->neighbors[i] == v) {
                memmove(g->neighbors + i,
                        g->neighbors + i + 1,
                        (size_t)(total - i - 1) * sizeof(int));
                break;
            }
        }
    }
    /* Shift all subsequent offsets by -1. */
    {
        int j;
        for (j = u + 1; j <= g->n; ++j) g->offsets[j] = g->offsets[j] - 1;
    }
}

/* Add undirected edge (u, v). Assumes it does not already exist. */
void ws_add_edge(WSGraph *g, int u, int v) {
    ws_insert_neighbor(g, u, v);
    ws_insert_neighbor(g, v, u);
}

/* Remove undirected edge (u, v). Assumes it exists. */
void ws_remove_edge(WSGraph *g, int u, int v) {
    ws_erase_neighbor(g, u, v);
    ws_erase_neighbor(g, v, u);
}

/* Check if edge (u, v) exists using binary search on u's list. */
int ws_has_edge(WSGraph *g, int u, int v) {
    int start = g->offsets[u];
    int end   = g->offsets[u+1];
    return ws_bs_has(g->neighbors + start, end - start, v);
}

/* ============================================================
 * 2) PRNG (LCG) — simple and readable (no bitwise)
 * ============================================================ */

/*
 * A simple linear-congruential generator (LCG).
 * Parameters are common textbook values; this is not cryptographically secure.
 * State must be non-zero. If the user seeds with 0, we promote it to 1.
 *
 * next = (a * current + c) mod m
 * We pick: a = 1103515245, c = 12345, m = 2^31  (implemented via % 2147483648)
 */
void ws_prng_seed(WSPrng *rng, int seed) {
    if (rng == 0) return;
    if (seed == 0) seed = 1;
    rng->state = seed;
}

int ws_prng_next_int_raw(WSPrng *rng) {
    /* Keep state in [0, 2^31) to avoid negatives; use long for overflow safety. */
    long a = 1103515245L;
    long c = 12345L;
    long m = 2147483648L; /* 2^31 */

    long s = (long)rng->state;
    s = (a * s + c) % m;
    rng->state = (int)s;
    return rng->state;
}

/* Float in [0,1). */
float ws_prng_next(WSPrng *rng) {
    if (rng == 0) return 0.0f;
    /* Divide by 2^31 to scale into [0,1). */
    int r = ws_prng_next_int_raw(rng);
    return (float)r / 2147483648.0f;
}

/* Integer in [0, max). If max == 0, return 0. */
int ws_prng_nextInt(WSPrng *rng, int max) {
    if (rng == 0 || max <= 0) return 0;
    /* Use float scaling; simple and portable. */
    float f = ws_prng_next(rng);
    int v = (int)(f * (float)max);
    if (v >= max) v = max - 1; /* guard against rounding edge */
    return v;
}

/* ============================================================
 * 3) Public: Build ring lattice
 * ============================================================ */

int ws_build_ring(int n, int k, WSGraph *out) {
    int u, d;

    if (out == 0) return -10;
    if (n <= 0 || (k % 2) != 0 || k >= n) return -1;

    /* Start with a clean struct. */
    out->n = n;
    out->k = k;
    out->offsets  = 0;
    out->neighbors= 0;
    out->capacity = 0;

    /* Allocate offsets (n+1) and neighbors (n*k). */
    out->offsets = (int*)calloc((size_t)(n + 1), sizeof(int));
    if (out->offsets == 0) return -2;

    out->capacity = n * k; /* total neighbor entries = n*k */
    out->neighbors = (int*)malloc((size_t)out->capacity * sizeof(int));
    if (out->neighbors == 0) {
        free(out->offsets); out->offsets = 0;
        return -3;
    }

    /* Each node has degree k, so offsets[u] = u*k. */
    for (u = 0; u < n; ++u) out->offsets[u] = u * k;
    out->offsets[n] = n * k;

    /* Fill each node's neighbor list (sorted), k/2 on each side around the ring. */
    {
        int half = k / 2;
        for (u = 0; u < n; ++u) {
            /* Collect neighbors into a small temporary array. */
            int tmp_count = 0;
            /* If you expect k larger than ~128, enlarge this buffer. */
            int tmp[128];

            for (d = 1; d <= half; ++d) {
                int a = (u + d) % n;
                int b = (u + n - d) % n;
                tmp[tmp_count++] = a;
                tmp[tmp_count++] = b;
            }

            /* Insertion sort is fine because k is typically small. */
            {
                int i;
                for (i = 1; i < tmp_count; ++i) {
                    int x = tmp[i];
                    int j = i;
                    while (j > 0 && tmp[j - 1] > x) {
                        tmp[j] = tmp[j - 1];
                        j = j - 1;
                    }
                    tmp[j] = x;
                }
            }

            /* Copy sorted neighbors into the main array region for node u. */
            memcpy(out->neighbors + out->offsets[u],
                   tmp,
                   (size_t)k * sizeof(int));
        }
    }

    return 0;
}

/* ============================================================
 * 4) Public: Rewiring pass
 * ============================================================ */

int ws_rewire_pass(WSGraph *g, WSPrng *rng, float beta, float frac) {
    int E, i, changed = 0;

    if (g == 0 || rng == 0) return -1;
    if (beta <= 0.0f || frac <= 0.0f) return 0;

    E = ws_num_edges(g);
    if (E <= 0) return 0;

    /* Target rewires ~ frac * beta * E (rounded stochastically). */
    {
        float expected = frac * beta * (float)E;
        int to_do = (int)expected;
        float rem = expected - (float)to_do;
        if (rem > 0.0f && ws_prng_next(rng) < rem) to_do = to_do + 1;

        if (to_do <= 0) return 0;

        /* Snapshot current unique edges as pairs (u < v) so we can safely mutate. */
        {
            int *elist = (int*)malloc((size_t)(2 * E) * sizeof(int));
            int m = 0;
            int u;

            if (elist == 0) return -2;

            for (u = 0; u < g->n; ++u) {
                int s = g->offsets[u];
                int e = g->offsets[u + 1];
                for (i = s; i < e; ++i) {
                    int v = g->neighbors[i];
                    if (u < v) {
                        elist[2 * m + 0] = u;
                        elist[2 * m + 1] = v;
                        m = m + 1;
                    }
                }
            }

            /* Apply to_do rewires. Keep it simple: try random new endpoints. */
            for (i = 0; i < to_do; ++i) {
                if (m <= 0) break;

                /* Pick a random existing edge (u, v) with u < v. */
                {
                    int pick = ws_prng_nextInt(rng, m);
                    int u0 = elist[2 * pick + 0];
                    int v0 = elist[2 * pick + 1];

                    /* Attempt to rewire u0's end to a new w. */
                    int attempt = 0;
                    int applied = 0;

                    while (attempt < 100) {
                        int w = ws_prng_nextInt(rng, g->n);
                        attempt = attempt + 1;

                        if (w == u0 || w == v0) continue;    /* avoid self loop */
                        if (!ws_has_edge(g, u0, v0)) break;   /* edge may be gone */

                        if (ws_has_edge(g, u0, w)) continue;  /* avoid duplicate */

                        /* Apply: remove (u0, v0), add (u0, w). */
                        ws_remove_edge(g, u0, v0);
                        ws_add_edge(g, u0, w);
                        changed = changed + 1;
                        applied = 1;
                        break;
                    }
                    (void)applied; /* kept for readability; unused value ok */
                }
            }

            free(elist);
        }
    }

    return changed;
}

/* ============================================================
 * 5) Public: Export edge list (unique pairs)
 * ============================================================ */

int ws_export_edgelist(WSGraph *g, int *out_pairs, int max_pairs) {
    int u, i, w = 0;
    int E;

    if (g == 0 || out_pairs == 0) return -1;

    E = ws_num_edges(g);
    if (max_pairs < E) return -2;

    for (u = 0; u < g->n; ++u) {
        int s = g->offsets[u];
        int e = g->offsets[u + 1];
        for (i = s; i < e; ++i) {
            int v = g->neighbors[i];
            if (u < v) {
                out_pairs[2 * w + 0] = u;
                out_pairs[2 * w + 1] = v;
                w = w + 1;
            }
        }
    }
    return w;  /* number of edges */
}

/* ============================================================
 * 6) Public: Free graph storage
 * ============================================================ */

void ws_free(WSGraph *g) {
    if (g == 0) return;
    free(g->offsets);   g->offsets   = 0;
    free(g->neighbors); g->neighbors = 0;
    g->n = 0; g->k = 0; g->capacity = 0;
}

/* ============================================================
 * 7) Public: Clustering coefficient computation
 * ============================================================ */

/*
 * Count triangles at node u using two-pointer set intersection.
 *
 * For each neighbor 'a' of u, count common neighbors between:
 *   - N(a): neighbors of a
 *   - N(u) in positions after 'a' (to avoid double counting)
 *
 * Since neighbor arrays are sorted, we use two-pointer technique
 * for O(deg^2) complexity per node.
 */
int ws_count_triangles_at_node(WSGraph *g, int u) {
    int start_u, end_u, deg_u;
    int i, triangles;
    int *Nu;

    if (g == 0) return 0;

    start_u = g->offsets[u];
    end_u = g->offsets[u + 1];
    deg_u = end_u - start_u;

    if (deg_u < 2) return 0;

    Nu = g->neighbors + start_u;
    triangles = 0;

    /* For each neighbor 'a' of u */
    for (i = 0; i < deg_u; ++i) {
        int a = Nu[i];
        int start_a = g->offsets[a];
        int end_a = g->offsets[a + 1];
        int deg_a = end_a - start_a;
        int *Na = g->neighbors + start_a;

        /* Two-pointer intersection: Na vs Nu[i+1..deg_u-1] */
        int p1 = 0;           /* pointer in Na */
        int p2 = i + 1;       /* pointer in Nu (only neighbors after index i) */

        while (p1 < deg_a && p2 < deg_u) {
            int node_a = Na[p1];
            int node_b = Nu[p2];

            if (node_a == node_b) {
                /* Found common neighbor -> triangle (u, a, common) */
                triangles = triangles + 1;
                p1 = p1 + 1;
                p2 = p2 + 1;
            } else if (node_a < node_b) {
                p1 = p1 + 1;
            } else {
                p2 = p2 + 1;
            }
        }
    }

    return triangles;
}

/*
 * Compute average clustering coefficient (Watts-Strogatz formula).
 *
 * For each node i with degree k_i >= 2:
 *   C_i = 2 * triangles_i / (k_i * (k_i - 1))
 *
 * Returns: C = (1/N) * sum(C_i) for nodes with degree >= 2.
 */
double ws_avg_clustering(WSGraph *g) {
    int u;
    double sum_clustering;
    int count_nodes;

    if (g == 0) return 0.0;
    if (g->n <= 0) return 0.0;

    sum_clustering = 0.0;
    count_nodes = 0;

    /* Compute C_i for each node with degree >= 2 */
    for (u = 0; u < g->n; ++u) {
        int start = g->offsets[u];
        int end = g->offsets[u + 1];
        int degree = end - start;

        if (degree >= 2) {
            int triangles = ws_count_triangles_at_node(g, u);
            int denominator = degree * (degree - 1);
            double Ci = (2.0 * (double)triangles) / (double)denominator;

            sum_clustering = sum_clustering + Ci;
            count_nodes = count_nodes + 1;
        }
        /* Nodes with degree < 2 have C_i = 0 (not included in average) */
    }

    /* Average: (1/N) * sum(C_i) */
    if (count_nodes > 0) {
        return sum_clustering / (double)count_nodes;
    }

    return 0.0;
}
