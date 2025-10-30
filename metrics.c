/* metrics.c  —  textbook, readable implementations of core graph metrics.
 *
 * Computes:
 *  - Connected components (IDs, sizes, giant, reachable pairs)
 *  - Average local clustering (Watts–Strogatz definition)
 *  - Global clustering / transitivity (optional helper)
 *  - Exact Average Path Length on the giant component
 *  - Effective diameter D90 (90th percentile distance on GC)
 *  - Small-world index sigma from precomputed baselines
 *
 * Assumptions:
 *  - Graph is undirected and simple (no multi-edges, no self-loops).
 *  - graph_neighbors() returns each vertex’s neighbor list in ASCENDING order.
 */

#include "metrics.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* --------------------------- Components (BFS) --------------------------- */

void components_init(Components* comp, int n) {
    if (!comp) return;
    comp->num_components = 0;
    comp->comp_id    = (int*)  malloc((size_t)n * sizeof(int));
    comp->comp_sizes = (int*)  calloc((size_t)n, sizeof(int));
    comp->giant_size = 0;
    comp->total_pairs = 0;
}

void components_free(Components* comp) {
    if (!comp) return;
    free(comp->comp_id);    comp->comp_id = NULL;
    free(comp->comp_sizes); comp->comp_sizes = NULL;
    comp->num_components = 0;
    comp->giant_size = 0;
    comp->total_pairs = 0;
}

/* Label connected components with a standard BFS; also compute:
 *  - comp_sizes[cid]
 *  - giant_size
 *  - total_pairs = sum over components of (size choose 2) = #reachable unordered pairs
 */
void find_components(Graph* g, Components* comp) {
    if (!g || !comp || g->n <= 0) return;

    char* visited = (char*)calloc((size_t)g->n, sizeof(char));
    int*  queue   = (int*) malloc((size_t)g->n * sizeof(int));

    comp->num_components = 0;
    comp->giant_size = 0;
    comp->total_pairs = 0;

    for (int s = 0; s < g->n; s++) {
        if (visited[s]) continue;

        int cid = comp->num_components++;
        int csize = 0;

        int qh = 0, qt = 0;
        queue[qt++] = s;
        visited[s] = 1;
        comp->comp_id[s] = cid;

        while (qh < qt) {
            int u = queue[qh++];
            csize++;

            int deg = 0;
            int* nbrs = graph_neighbors(g, u, &deg);
            for (int i = 0; i < deg; i++) {
                int v = nbrs[i];
                if (!visited[v]) {
                    visited[v] = 1;
                    comp->comp_id[v] = cid;
                    queue[qt++] = v;
                }
            }
        }

        comp->comp_sizes[cid] = csize;
        if (csize > comp->giant_size) comp->giant_size = csize;
        if (csize > 1) comp->total_pairs += (long long)csize * (csize - 1) / 2;
    }

    free(visited);
    free(queue);
}

/* -------------------- Clustering (local & global) ---------------------- */

/* Count e_i = # edges among neighbors of u.
 * Neighbor lists are sorted; we intersect N(u) with N(v) for v in N(u) and only count
 * pairs (v,w) with w > v to avoid double counting. This counts each neighbor pair once.
 */
static int count_neighbor_links(Graph* g, int u) {
    int deg_u = 0;
    int* Nu = graph_neighbors(g, u, &deg_u);
    if (deg_u < 2) return 0;

    int e_i = 0;

    for (int idx = 0; idx < deg_u; idx++) {
        int v = Nu[idx];

        int deg_v = 0;
        int* Nv = graph_neighbors(g, v, &deg_v);

        /* Two-pointer intersection between Nv and the suffix of Nu (elements after v).
           We only consider neighbors w in Nu with w > v to count each pair once. */
        int pi = idx + 1;   /* pointer in Nu */
        int pv = 0;         /* pointer in Nv */
        while (pi < deg_u && pv < deg_v) {
            int a = Nu[pi];
            int b = Nv[pv];
            if (a == b) { e_i++; pi++; pv++; }
            else if (a < b) { pi++; }
            else { pv++; }
        }
    }
    return e_i;  /* unordered neighbor–neighbor edges */
}

/* Average LOCAL clustering: C = (1/N) * sum_i Ci
 * with Ci = 2*e_i / (k_i*(k_i-1)), and Ci=0 if k_i<2.
 */
double compute_clustering_local(Graph* g) {
    if (!g || g->n <= 0) return 0.0;

    double sumCi = 0.0;

    for (int u = 0; u < g->n; u++) {
        int k = graph_degree(g, u);
        if (k < 2) { /* Ci = 0 */ continue; }

        int ei = count_neighbor_links(g, u);     /* edges among neighbors of u */
        int denom = k * (k - 1) / 2;             /* # possible neighbor–neighbor edges */
        if (denom > 0) {
            double Ci = (double)ei / (double)denom;
            sumCi += Ci;
        }
    }
    return sumCi / (double)g->n;  /* nodes with k<2 contribute 0 implicitly */
}

/* Global clustering / transitivity:
 *   C_global = (3 * #triangles) / #wedges
 * where #wedges = sum_i k_i*(k_i-1)/2 and each triangle is counted once overall.
 * We compute triangles via per-node e_i and then sum: sum_i e_i = 3 * #triangles.
 */
double compute_clustering_global(Graph* g) {
    if (!g || g->n <= 0) return 0.0;

    long long sum_ei = 0;      /* sum over nodes of e_i */
    long long wedges = 0;      /* sum over nodes of k_i*(k_i-1)/2 */

    for (int u = 0; u < g->n; u++) {
        int k = graph_degree(g, u);
        if (k >= 2) {
            sum_ei += (long long)count_neighbor_links(g, u);
            wedges += (long long)k * (k - 1) / 2;
        }
    }
    if (wedges == 0) return 0.0;
    /* 3*triangles / wedges  == (3*(sum_ei/3)) / wedges == sum_ei / wedges */
    return (double)sum_ei / (double)wedges;
}

/* ----------------- Exact APL + Effective Diameter (GC) ------------------ */

/* BFS from a single source over the giant component only.
 * Fills dist[] with hop counts (−1 = unreachable).
 */
static void bfs_from_source(Graph* g, int src, int* dist) {
    int n = g->n;
    int* q = (int*)malloc((size_t)n * sizeof(int));
    int qh = 0, qt = 0;

    for (int i = 0; i < n; i++) dist[i] = -1;
    dist[src] = 0;
    q[qt++] = src;

    while (qh < qt) {
        int u = q[qh++];

        int deg = 0;
        int* nbrs = graph_neighbors(g, u, &deg);
        for (int i = 0; i < deg; i++) {
            int v = nbrs[i];
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                q[qt++] = v;
            }
        }
    }
    free(q);
}

/* Compute:
 *  - avg_path_len (ADS/APL) over all reachable unordered pairs in the GIANT COMPONENT
 *  - eff_d90       the 90th percentile distance (effective diameter)
 *  - reachable_pairs (optional out) = number of unordered pairs actually averaged
 *
 * Requires a Components labelling already computed; we detect the GC id first.
 */
void compute_apl_and_d90_on_gc(Graph* g,
                               const Components* comp,
                               double* avg_path_len,
                               double* eff_d90,
                               long long* reachable_pairs)
{
    if (!g || g->n <= 0 || !comp || comp->num_components <= 0) {
        if (avg_path_len) *avg_path_len = 0.0;
        if (eff_d90) *eff_d90 = 0.0;
        if (reachable_pairs) *reachable_pairs = 0;
        return;
    }

    /* Identify the giant component ID (first max). */
    int gc_id = -1, gc_size = 0;
    for (int cid = 0; cid < comp->num_components; cid++) {
        if (comp->comp_sizes[cid] > gc_size) {
            gc_size = comp->comp_sizes[cid];
            gc_id = cid;
        }
    }
    if (gc_id < 0 || gc_size < 2) {
        if (avg_path_len) *avg_path_len = 0.0;
        if (eff_d90) *eff_d90 = 0.0;
        if (reachable_pairs) *reachable_pairs = 0;
        return;
    }

    /* Prepare a list of nodes in the giant component for iteration. */
    int* gc_nodes = (int*)malloc((size_t)gc_size * sizeof(int));
    int gi = 0;
    for (int u = 0; u < g->n; u++) if (comp->comp_id[u] == gc_id) gc_nodes[gi++] = u;

    /* Distance histogram up to (gc_size-1). */
    long long* hist = (long long*)calloc((size_t)gc_size, sizeof(long long));
    long long pair_count = 0;
    long long dist_sum   = 0;

    int* dist = (int*)malloc((size_t)g->n * sizeof(int));

    /* BFS from every GC node; accumulate distances to nodes with index > sourceIndex
       to count each unordered pair exactly once and avoid double counting. */
    for (int idx = 0; idx < gc_size; idx++) {
        int s = gc_nodes[idx];
        bfs_from_source(g, s, dist);

        for (int j = idx + 1; j < gc_size; j++) {
            int t = gc_nodes[j];
            int d = dist[t];
            if (d > 0) { /* (s,t) is reachable and not the same node */
                dist_sum += d;
                pair_count++;
                if (d < gc_size) hist[d]++;     /* d is <= (gc_size-1) */
            }
        }
    }

    /* Average path length on GC */
    double L = (pair_count > 0) ? (double)dist_sum / (double)pair_count : 0.0;

    /* Effective diameter D90: smallest d with cumulative ≥ 0.9 * pair_count */
    double D90 = 0.0;
    if (pair_count > 0) {
        long long target = (long long)ceil(0.90 * (double)pair_count);
        long long cum = 0;
        for (int d = 1; d < gc_size; d++) {
            cum += hist[d];
            if (cum >= target) { D90 = (double)d; break; }
        }
    }

    if (avg_path_len)   *avg_path_len = L;
    if (eff_d90)        *eff_d90 = D90;
    if (reachable_pairs) *reachable_pairs = pair_count;

    free(dist);
    free(hist);
    free(gc_nodes);
}

/* ---------------------- Small-world index helper ----------------------- */

double small_world_index_sigma(double C, double L, double C0, double L0) {
    /* Guard against zero baselines */
    if (C0 <= 0.0 || L0 <= 0.0) return 0.0;
    double num = (C / C0);
    double den = (L / L0);
    if (den == 0.0) return 0.0;
    return num / den;
}

/* Wrapper for backward compatibility and to match the API expected by sim.c
 * By default, uses the local (Watts-Strogatz) definition.
 * To use global transitivity instead, change this to call compute_clustering_global(g).
 */
double compute_clustering(Graph* g) {
    return compute_clustering_local(g);  /* Default: WS local definition */
}
