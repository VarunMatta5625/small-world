#include "sim.h"
#include "graph.h"
#include "rng.h"
#include "layout.h"
#include "bfs.h"
#include "metrics.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

/* Limits for visualization buffers */
#define MAX_PATH_EDGES   4096
#define MAX_CHART_POINTS 10000

/* A single time-series datapoint for charts. */
typedef struct {
    int tick;
    double ads;
    double clustering;
    double giant_share;
    double compute_ms;
} ChartPoint;

/* Global simulation state (intentionally not static per your constraints). */
struct SimState {
    Graph graph;
    Layout layout;
    BFS bfs;
    RNG rng;
    Components components;

    /* Parameters */
    int n;
    int k;
    float beta;
    int seed;
    int layout_type; /* 0=ring, 1=force */

    /* ADS running stats */
    long long pairs_seen;
    long long dist_sum;
    double last_ads;

    /* Pair traversal state */
    int current_src;
    int current_dst;
    int current_dist;
    int bfs_done;       /* BFS from current_src started/completed? */
    int path_ready;     /* Is a path prepared for (src,dst)? */

    int path_edges[MAX_PATH_EDGES * 2]; /* flattened (u,v) pairs */
    int path_edge_count;
    int path_edges_revealed;

    /* All edges for drawing */
    int* all_edges;     /* flattened (u,v) pairs */
    int  all_edges_count;

    /* Metrics */
    double clustering;
    double giant_share;

    /* Timing */
    double compute_ms;

    /* Chart data */
    ChartPoint* chart_points;
    int chart_count;
    int chart_capacity;

    /* Tick counter */
    int tick_count;

    /* Viewport (kept inside layout as well) */
} sim;

/* Forward declarations for helper functions (not static per request). */
void sim_prepare_path(void);
void sim_advance_pair(void);

EMSCRIPTEN_KEEPALIVE void sim_init(int width, int height) {
    memset(&sim, 0, sizeof(sim));

    /* Defaults */
    sim.n = 140;
    sim.k = 6;
    sim.beta = 0.2f;
    sim.seed = 12345;
    sim.layout_type = 0; /* ring */

    graph_init(&sim.graph, sim.n);
    layout_init(&sim.layout, sim.n, LAYOUT_RING, (float)width, (float)height);
    bfs_init(&sim.bfs, &sim.graph);
    rng_init(&sim.rng, sim.seed);
    components_init(&sim.components, sim.n);

    sim.chart_capacity = MAX_CHART_POINTS;
    sim.chart_points = (ChartPoint*)malloc((size_t)sim.chart_capacity * sizeof(ChartPoint));
    sim.chart_count = 0;

    sim.all_edges = 0;
    sim.all_edges_count = 0;

    sim.tick_count = 0;
}

EMSCRIPTEN_KEEPALIVE void sim_free(void) {
    graph_free(&sim.graph);
    layout_free(&sim.layout);
    bfs_free(&sim.bfs);
    components_free(&sim.components);
    if (sim.chart_points) free(sim.chart_points);
    if (sim.all_edges) free(sim.all_edges);
    memset(&sim, 0, sizeof(sim));
}

EMSCRIPTEN_KEEPALIVE void sim_set_params(
    int n,
    int k,
    float beta01,
    int seed,
    int layout_type
) {
    sim.n = n;
    sim.k = k;
    sim.beta = beta01;
    sim.seed = seed;
    sim.layout_type = layout_type;
}

EMSCRIPTEN_KEEPALIVE void sim_build_graph(void) {
    /* Free old and reinit with current sizes */
    graph_free(&sim.graph);
    layout_free(&sim.layout);
    bfs_free(&sim.bfs);
    components_free(&sim.components);

    graph_init(&sim.graph, sim.n);
    layout_init(&sim.layout, sim.n,
                (sim.layout_type == 0 ? LAYOUT_RING : LAYOUT_FORCE),
                sim.layout.width, sim.layout.height);
    bfs_init(&sim.bfs, &sim.graph);
    rng_init(&sim.rng, sim.seed);
    components_init(&sim.components, sim.n);

    /* Build new WS graph */
    ws_build(&sim.graph, sim.n, sim.k, sim.beta, &sim.rng);

    /* Build edge list for drawing (unique u<v) */
    if (sim.all_edges) free(sim.all_edges);
    sim.all_edges_count = sim.graph.m;
    sim.all_edges = (int*)malloc((size_t)sim.all_edges_count * 2 * sizeof(int));

    int edge_idx = 0;
    for (int u = 0; u < sim.graph.n; u++) {
        int deg = 0;
        int* nbrs = graph_neighbors(&sim.graph, u, &deg);
        for (int i = 0; i < deg; i++) {
            int v = nbrs[i];
            if (u < v) {
                sim.all_edges[edge_idx * 2 + 0] = u;
                sim.all_edges[edge_idx * 2 + 1] = v;
                edge_idx++;
            }
        }
    }

    /* Layout */
    if (sim.layout_type == 0) layout_compute_ring(&sim.layout);
    else layout_init_force(&sim.layout, &sim.rng);

    /* Components and metrics */
    find_components(&sim.graph, &sim.components);
    sim.clustering = compute_clustering(&sim.graph);
    sim.giant_share = (sim.n > 0) ? ((double)sim.components.giant_size / (double)sim.n) : 0.0;

    /* Reset ADS state */
    sim.pairs_seen = 0;
    sim.dist_sum = 0;
    sim.last_ads = 0.0;

    sim.current_src = 0;
    sim.current_dst = 0;
    sim.current_dist = 0;
    sim.bfs_done = 0;
    sim.path_ready = 0;
    sim.path_edge_count = 0;
    sim.path_edges_revealed = 0;

    sim.tick_count = 0;
    sim.chart_count = 0;
}

/* Prepare path display for current (src,dst). */
void sim_prepare_path(void) {
    if (sim.current_dst >= sim.n) {
        sim.path_ready = 0;
        sim.path_edge_count = 0;
        return;
    }

    /* Unreachable? In BFS we used BFS_INF=2147483647, so check against that. */
    if (sim.bfs.dist[sim.current_dst] >= 2147483647) {
        sim.path_ready = 0;
        sim.path_edge_count = 0;
        sim.current_dist = 2147483647;
        return;
    }

    sim.current_dist = sim.bfs.dist[sim.current_dst];
    sim.path_edge_count = bfs_path(&sim.bfs, sim.current_src, sim.current_dst, sim.path_edges);
    sim.path_edges_revealed = 0;
    sim.path_ready = 1;
}

/* Advance to next (src,dst) pair in lexicographic order. */
void sim_advance_pair(void) {
    sim.current_dst++;
    if (sim.current_dst >= sim.n) {
        sim.current_src++;
        if (sim.current_src >= sim.n) {
            /* Done over all pairs */
            sim.current_src = sim.n;
            sim.current_dst = sim.n;
            sim.bfs_done = 1;
            sim.path_ready = 0;
            return;
        }
        bfs_run(&sim.bfs, sim.current_src);
        sim.bfs_done = 1;
        sim.current_dst = sim.current_src + 1;
    }
    if (sim.current_dst < sim.n) sim_prepare_path();
}

EMSCRIPTEN_KEEPALIVE void sim_tick(int layout_steps, int path_edges_per_tick) {
    clock_t start = clock();

    /* Layout */
    if (sim.layout_type == 1) {
        for (int i = 0; i < layout_steps; i++) layout_step_force(&sim.layout, &sim.graph);
    }

    /* Pair-wise traversal visualization */
    if (sim.current_src >= sim.n) {
        /* all pairs processed */
        clock_t end2 = clock();
        sim.compute_ms = ((double)(end2 - start) / CLOCKS_PER_SEC) * 1000.0;
        sim.tick_count++;
        return;
    }

    /* Start first BFS if needed */
    if (!sim.bfs_done && sim.current_src == 0 && sim.current_dst == 0) {
        bfs_run(&sim.bfs, sim.current_src);
        sim.bfs_done = 1;
        sim.current_dst = sim.current_src + 1;
        if (sim.current_dst < sim.n) sim_prepare_path();
    }

    /* Reveal a few edges of the current path */
    if (sim.path_ready && sim.path_edges_revealed < sim.path_edge_count) {
        int to_rev = path_edges_per_tick;
        if (sim.path_edges_revealed + to_rev > sim.path_edge_count) {
            to_rev = sim.path_edge_count - sim.path_edges_revealed;
        }
        sim.path_edges_revealed += to_rev;

        /* If fully revealed, account the distance and move on */
        if (sim.path_edges_revealed >= sim.path_edge_count) {
            if (sim.current_dist < 2147483647) {
                sim.pairs_seen += 1;
                sim.dist_sum += (long long)sim.current_dist;
                if (sim.pairs_seen > 0) {
                    sim.last_ads = (double)sim.dist_sum / (double)sim.pairs_seen;
                }

                if (sim.chart_count < MAX_CHART_POINTS) {
                    sim.chart_points[sim.chart_count].tick = sim.tick_count;
                    sim.chart_points[sim.chart_count].ads = sim.last_ads;
                    sim.chart_points[sim.chart_count].clustering = sim.clustering;
                    sim.chart_points[sim.chart_count].giant_share = sim.giant_share;
                    sim.chart_points[sim.chart_count].compute_ms = sim.compute_ms;
                    sim.chart_count++;
                }
            }
            sim_advance_pair();
        }
    } else if (!sim.path_ready) {
        /* If no path to reveal (unreachable or zero-length), advance immediately. */
        sim_advance_pair();
    }

    clock_t end = clock();
    sim.compute_ms = ((double)(end - start) / CLOCKS_PER_SEC) * 1000.0;
    sim.tick_count++;
}

/* Accessors for drawing and stats */
EMSCRIPTEN_KEEPALIVE int sim_n(void) { return sim.graph.n; }
EMSCRIPTEN_KEEPALIVE int sim_m(void) { return sim.graph.m; }
EMSCRIPTEN_KEEPALIVE float* sim_pos_x(void) { return sim.layout.pos_x; }
EMSCRIPTEN_KEEPALIVE float* sim_pos_y(void) { return sim.layout.pos_y; }

EMSCRIPTEN_KEEPALIVE int   sim_path_edge_count(void) { return sim.path_edges_revealed; }
EMSCRIPTEN_KEEPALIVE int*  sim_path_edges_ptr(void)  { return sim.path_edges; }

EMSCRIPTEN_KEEPALIVE int sim_current_src(void)  { return sim.current_src; }
EMSCRIPTEN_KEEPALIVE int sim_current_dst(void)  { return sim.current_dst; }
EMSCRIPTEN_KEEPALIVE int sim_current_dist(void) { return sim.current_dist; }

EMSCRIPTEN_KEEPALIVE double sim_ads(void)         { return sim.last_ads; }
EMSCRIPTEN_KEEPALIVE double sim_clustering(void)  { return sim.clustering; }
EMSCRIPTEN_KEEPALIVE double sim_giant_share(void) { return sim.giant_share; }
EMSCRIPTEN_KEEPALIVE double sim_compute_ms(void)  { return sim.compute_ms; }

EMSCRIPTEN_KEEPALIVE double sim_progress(void) {
    if (sim.components.total_pairs <= 0) return 100.0;
    return (double)sim.pairs_seen / (double)sim.components.total_pairs * 100.0;
}

EMSCRIPTEN_KEEPALIVE int   sim_all_edges_count(void) { return sim.all_edges_count; }
EMSCRIPTEN_KEEPALIVE int*  sim_all_edges_ptr(void)   { return sim.all_edges; }

/* CSV export of time-series chart points. */
EMSCRIPTEN_KEEPALIVE int sim_export_csv(char* out_buf, int cap) {
    if (!out_buf || cap <= 0) return 0;

    const char* header = "tick,ads,clustering,giant,compute_ms\n";
    int pos = 0;

    /* Write header */
    int wrote = snprintf(out_buf + pos, (size_t)(cap - pos), "%s", header);
    if (wrote < 0) return 0;
    pos += wrote;

    /* Write rows */
    for (int i = 0; i < sim.chart_count && pos < cap - 128; i++) {
        wrote = snprintf(out_buf + pos, (size_t)(cap - pos),
                         "%d,%.6f,%.6f,%.6f,%.3f\n",
                         sim.chart_points[i].tick,
                         sim.chart_points[i].ads,
                         sim.chart_points[i].clustering,
                         sim.chart_points[i].giant_share,
                         sim.chart_points[i].compute_ms);
        if (wrote < 0) break;
        pos += wrote;
    }
    return pos;
}
