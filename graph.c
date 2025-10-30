#include "graph.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================
 * Utility helpers (no static, no bitwise)
 * ============================================================ */

/* Binary search in sorted array a[0..len-1]; returns 1 if found, else 0. */
int g_bin_has(int* a, int len, int x) {
    int lo = 0, hi = len;
    while (lo < hi) {
        int m = lo + (hi - lo) / 2;
        int v = a[m];
        if (v == x) return 1;
        if (v < x)  lo = m + 1; else hi = m;
    }
    return 0;
}

/* Insert x into sorted region a[start..end-1] by shifting tail. */
void g_sorted_insert(Graph* g, int u, int x) {
    int start = g->offs[u];
    int end   = g->offs[u + 1];
    int len   = end - start;

    /* find insert position */
    int lo = 0, hi = len;
    while (lo < hi) {
        int m = lo + (hi - lo) / 2;
        if (g->nbrs[start + m] < x) lo = m + 1; else hi = m;
    }
    int ins_pos = start + lo;

    /* total used slots */
    int total = g->offs[g->n];

    /* ensure capacity */
    /* grow neighbors by doubling */
    if (!g->nbrs) {
        g->nbrs = (int*)malloc(16 * sizeof(int));
    }
    /* compute current capacity from m (2*m) or store separately. For simplicity:
       realloc to total+1 elements each time we need (still readable, small). */
    g->nbrs = (int*)realloc(g->nbrs, (size_t)(total + 1) * sizeof(int));

    /* move suffix */
    memmove(g->nbrs + ins_pos + 1,
            g->nbrs + ins_pos,
            (size_t)(total - ins_pos) * sizeof(int));

    /* insert */
    g->nbrs[ins_pos] = x;

    /* shift tail offsets */
    for (int i = u + 1; i <= g->n; i++) g->offs[i] = g->offs[i] + 1;
}

/* Erase x from u's neighbor list (assumes exists). */
void g_sorted_erase(Graph* g, int u, int x) {
    int start = g->offs[u];
    int end   = g->offs[u + 1];
    int total = g->offs[g->n];

    for (int i = start; i < end; i++) {
        if (g->nbrs[i] == x) {
            /* shift left the suffix */
            memmove(g->nbrs + i,
                    g->nbrs + i + 1,
                    (size_t)(total - i - 1) * sizeof(int));
            break;
        }
    }

    /* reduce offsets for tails */
    for (int j = u + 1; j <= g->n; j++) g->offs[j] = g->offs[j] - 1;

    /* shrink allocation a little is optional; skip for simplicity */
}

/* ============================================================
 * Basic API
 * ============================================================ */

void graph_init(Graph* g, int n) {
    if (!g) return;
    g->n = n;
    g->m = 0;
    g->offs = (int*)calloc((size_t)(n + 1), sizeof(int));
    g->nbrs = 0;
}

void graph_free(Graph* g) {
    if (!g) return;
    if (g->offs) free(g->offs);
    if (g->nbrs) free(g->nbrs);
    g->offs = 0;
    g->nbrs = 0;
    g->n = 0;
    g->m = 0;
}

void graph_clear(Graph* g) {
    if (!g) return;
    g->m = 0;
    if (g->offs) memset(g->offs, 0, (size_t)(g->n + 1) * sizeof(int));
    if (g->nbrs) { free(g->nbrs); g->nbrs = 0; }
}

int graph_degree(Graph* g, int u) {
    if (!g || !g->offs || u < 0 || u >= g->n) return 0;
    return g->offs[u + 1] - g->offs[u];
}

int* graph_neighbors(Graph* g, int u, int* deg) {
    if (!g || !g->offs || !g->nbrs || u < 0 || u >= g->n) {
        if (deg) *deg = 0;
        return 0;
    }
    if (deg) *deg = g->offs[u + 1] - g->offs[u];
    return g->nbrs + g->offs[u];
}

int graph_has_edge(Graph* g, int u, int v) {
    if (!g || u < 0 || v < 0 || u >= g->n || v >= g->n || u == v) return 0;
    int deg = 0;
    int* nbrs = graph_neighbors(g, u, &deg);
    if (!nbrs) return 0;
    return g_bin_has(nbrs, deg, v);
}

/* Add/remove undirected edges while keeping lists sorted. */
void g_add_edge(Graph* g, int u, int v) {
    g_sorted_insert(g, u, v);
    g_sorted_insert(g, v, u);
}
void g_remove_edge(Graph* g, int u, int v) {
    g_sorted_erase(g, u, v);
    g_sorted_erase(g, v, u);
}

/* ============================================================
 * Watts–Strogatz build (ring + rewiring)
 * ============================================================ */
void ws_build(Graph* g, int n, int k, float beta, RNG* rng) {
    if (!g) return;

    /* Normalize parameters. */
    if (n <= 0) { graph_clear(g); return; }
    if (k % 2 != 0) k = k - 1;
    if (k < 0) k = 0;
    if (k >= n) k = (n > 2) ? 2 : 0;
    if (beta < 0.0f) beta = 0.0f;
    if (beta > 1.0f) beta = 1.0f;

    graph_clear(g);
    g->n = n;

    /* Preallocate ring lattice adjacency:
       Each node has degree k, so total neighbor slots = n*k. */
    if (g->offs) free(g->offs);
    g->offs = (int*)calloc((size_t)(n + 1), sizeof(int));
    if (g->nbrs) { free(g->nbrs); g->nbrs = 0; }

    for (int u = 0; u < n; u++) g->offs[u] = u * k;
    g->offs[n] = n * k;
    g->nbrs = (int*)malloc((size_t)(n * k) * sizeof(int));

    /* Fill sorted k-neighbor lists for ring. */
    int half = k / 2;
    for (int u = 0; u < n; u++) {
        /* Collect neighbors in small temporary array and insertion-sort it. */
        int tmp_count = 0;
        int tmp_cap = (k > 0 ? k : 1);
        int tmpbuf[256]; /* if k can exceed 256, switch to malloc; kept simple here. */

        /* neighbors on both sides around the ring */
        for (int d = 1; d <= half; d++) {
            int a = (u + d) % n;
            int b = (u + n - d) % n;
            tmpbuf[tmp_count++] = a;
            tmpbuf[tmp_count++] = b;
        }

        /* sort */
        for (int i = 1; i < tmp_count; i++) {
            int x = tmpbuf[i];
            int j = i;
            while (j > 0 && tmpbuf[j - 1] > x) {
                tmpbuf[j] = tmpbuf[j - 1];
                j = j - 1;
            }
            tmpbuf[j] = x;
        }

        /* copy to CSR range for u */
        if (k > 0) {
            memcpy(g->nbrs + g->offs[u], tmpbuf, (size_t)k * sizeof(int));
        }
    }

    g->m = (n * k) / 2; /* undirected edges */

    /* If beta is 0, we are done (pure ring). */
    if (beta <= 0.0f || k == 0) return;

    /* Rewiring: iterate approx beta fraction of edges, trying new endpoints. */
    int E = g->m;
    int target = (int)(beta * (float)E + 0.5f);

    /* Snapshot current unique edges as pairs (u < v). */
    int* edges = (int*)malloc((size_t)(2 * E) * sizeof(int));
    int m = 0;
    for (int u = 0; u < n; u++) {
        int deg = 0;
        int* nbrs = graph_neighbors(g, u, &deg);
        for (int i = 0; i < deg; i++) {
            int v = nbrs[i];
            if (u < v) { edges[2 * m + 0] = u; edges[2 * m + 1] = v; m++; }
        }
    }

    /* Randomly pick edges and try rewiring u->w. Keep symmetry and simplicity. */
    for (int t = 0; t < target && m > 0; t++) {
        int pick = rng_range(rng, m);
        int u = edges[2 * pick + 0];
        int v = edges[2 * pick + 1];

        /* try up to 100 candidates for w */
        int applied = 0;
        for (int attempt = 0; attempt < 100; attempt++) {
            int w = rng_range(rng, n);
            if (w == u || w == v) continue;
            if (!graph_has_edge(g, u, v)) break; /* edge might have been changed */
            if (graph_has_edge(g, u, w)) continue;

            g_remove_edge(g, u, v);
            g_add_edge(g, u, w);
            applied = 1;
            break;
        }
        (void)applied;
    }

    free(edges);
}
