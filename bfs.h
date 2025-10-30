#ifndef BFS_H
#define BFS_H
/*
 * Breadth-First Search (BFS) with predecessor tracking.
 * Plain C version:
 *   - uses int (no unsigned types)
 *   - no const/static
 *   - no bitwise operations
 *
 * The BFS object owns its working arrays and operates on a Graph
 * (declared in graph.h). Distances are counted in number of edges.
 */

#include "graph.h"

/* BFS working state. */
typedef struct {
    Graph* g;    /* graph to traverse */
    int* dist;   /* distance from source (INF = unreachable) */
    int* pred;   /* predecessor node to reconstruct shortest path tree */
    int* queue;  /* simple ring-less queue using two indices */
    char* visited; /* 0/1 flags for visited nodes */
} BFS;

/* Initialize BFS for graph g (allocates arrays of length g->n). */
void bfs_init(BFS* bfs, Graph* g);

/* Free internally allocated arrays. */
void bfs_free(BFS* bfs);

/* Run BFS from source node src. Fills dist[] and pred[] for reachable nodes. */
void bfs_run(BFS* bfs, int src);

/*
 * Reconstruct shortest path from src to dst using pred[] from last bfs_run.
 * Writes a forward list of edges (u,v) into path_out as pairs:
 *   path_out[0]=u0, path_out[1]=v0, path_out[2]=u1, path_out[3]=v1, ...
 * Returns number of edges in the path (0 if unreachable or src==dst).
 * Caller must provide a large enough buffer (up to n-1 edges).
 */
int bfs_path(BFS* bfs, int src, int dst, int* path_out);

#endif /* BFS_H */
