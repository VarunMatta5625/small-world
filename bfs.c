#include "bfs.h"
#include <stdlib.h>
#include <string.h>

/* Use large positive number to represent "unreachable". */
#define BFS_INF 2147483647

void bfs_init(BFS* bfs, Graph* g) {
    if (!bfs || !g) return;
    bfs->g = g;
    bfs->dist = (int*)malloc((size_t)g->n * sizeof(int));
    bfs->pred = (int*)malloc((size_t)g->n * sizeof(int));
    bfs->queue = (int*)malloc((size_t)g->n * sizeof(int));
    bfs->visited = (char*)calloc((size_t)g->n, sizeof(char));
}

void bfs_free(BFS* bfs) {
    if (!bfs) return;
    if (bfs->dist) free(bfs->dist);
    if (bfs->pred) free(bfs->pred);
    if (bfs->queue) free(bfs->queue);
    if (bfs->visited) free(bfs->visited);
    bfs->dist = 0;
    bfs->pred = 0;
    bfs->queue = 0;
    bfs->visited = 0;
    bfs->g = 0;
}

void bfs_run(BFS* bfs, int src) {
    if (!bfs || !bfs->g) return;

    Graph* g = bfs->g;
    int n = g->n;

    /* Initialize all arrays. */
    for (int i = 0; i < n; i++) {
        bfs->dist[i] = BFS_INF;
        bfs->pred[i] = -1;
        bfs->visited[i] = 0;
    }

    /* Queue managed by head (qh) and tail (qt). */
    int qh = 0, qt = 0;
    bfs->queue[qt++] = src;
    bfs->dist[src] = 0;
    bfs->pred[src] = src;
    bfs->visited[src] = 1;

    while (qh < qt) {
        int u = bfs->queue[qh++];

        int deg = 0;
        int* nbrs = graph_neighbors(g, u, &deg);

        for (int i = 0; i < deg; i++) {
            int v = nbrs[i];
            if (!bfs->visited[v]) {
                bfs->visited[v] = 1;
                bfs->dist[v] = bfs->dist[u] + 1;
                bfs->pred[v] = u;
                bfs->queue[qt++] = v;
            }
        }
    }
}

int bfs_path(BFS* bfs, int src, int dst, int* path_out) {
    if (!bfs || !bfs->g || !path_out) return 0;
    if (src == dst) return 0;
    if (dst < 0 || dst >= bfs->g->n) return 0;
    if (bfs->dist[dst] == BFS_INF) return 0;

    /* Reconstruct backwards into a temporary array of nodes. */
    int n = bfs->g->n;
    int* rnodes = (int*)malloc((size_t)n * sizeof(int));
    int len = 0;

    int cur = dst;
    while (cur != src && len < n) {
        rnodes[len++] = cur;
        cur = bfs->pred[cur];
        if (cur < 0) { free(rnodes); return 0; }
    }
    rnodes[len++] = src;

    /* Convert nodes to forward (u,v) edge pairs. */
    int edges = len - 1;
    for (int i = 0; i < edges; i++) {
        int u = rnodes[len - 1 - i];
        int v = rnodes[len - 2 - i];
        path_out[2 * i + 0] = u;
        path_out[2 * i + 1] = v;
    }

    free(rnodes);
    return edges;
}
