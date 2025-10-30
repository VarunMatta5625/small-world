#ifndef SIM_H
#define SIM_H
/*
 * Public WASM API for the Small-World (Watts–Strogatz) ADS simulation.
 * Plain C: int indices, no const/unsigned/static/bitwise.
 *
 * We keep EMSCRIPTEN_KEEPALIVE for exports. If building without Emscripten,
 * you can define EMSCRIPTEN_KEEPALIVE as empty in your build system.
 */

#include <emscripten.h>

/* Lifecycle */
EMSCRIPTEN_KEEPALIVE void sim_init(int width, int height);
EMSCRIPTEN_KEEPALIVE void sim_free(void);

/* Parameters & build */
EMSCRIPTEN_KEEPALIVE void sim_set_params(
    int n,
    int k,
    float beta01,    /* beta in [0,1] */
    int seed,
    int layout_type  /* 0=RING, 1=FORCE */
);
EMSCRIPTEN_KEEPALIVE void sim_build_graph(void);

/* Tick loop */
EMSCRIPTEN_KEEPALIVE void sim_tick(int layout_steps, int path_edges_per_tick);

/* Drawing data */
EMSCRIPTEN_KEEPALIVE int   sim_n(void);
EMSCRIPTEN_KEEPALIVE int   sim_m(void);
EMSCRIPTEN_KEEPALIVE float* sim_pos_x(void);
EMSCRIPTEN_KEEPALIVE float* sim_pos_y(void);

/* All graph edges for drawing (pairs u,v) */
EMSCRIPTEN_KEEPALIVE int   sim_all_edges_count(void);
EMSCRIPTEN_KEEPALIVE int*  sim_all_edges_ptr(void);

/* Path edges for current pair visualization */
EMSCRIPTEN_KEEPALIVE int   sim_path_edge_count(void);
EMSCRIPTEN_KEEPALIVE int*  sim_path_edges_ptr(void);

/* Current pair being processed */
EMSCRIPTEN_KEEPALIVE int sim_current_src(void);
EMSCRIPTEN_KEEPALIVE int sim_current_dst(void);
EMSCRIPTEN_KEEPALIVE int sim_current_dist(void);

/* Metrics */
EMSCRIPTEN_KEEPALIVE double sim_ads(void);
EMSCRIPTEN_KEEPALIVE double sim_clustering(void);
EMSCRIPTEN_KEEPALIVE double sim_giant_share(void);
EMSCRIPTEN_KEEPALIVE double sim_compute_ms(void);

/* Progress (0-100%) */
EMSCRIPTEN_KEEPALIVE double sim_progress(void);

/* Export CSV to out_buf (cap bytes). Returns bytes written. */
EMSCRIPTEN_KEEPALIVE int sim_export_csv(char* out_buf, int cap);

#endif /* SIM_H */
