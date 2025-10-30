// graph_metrics_omp.c
// Textbook-style graph metrics with OpenMP parallelism.
// Computes: |V|, |E|, avg degree, density, components (weak), average shortest path (reachable pairs),
// diameter, average local clustering, total triangles (undirected).
//
// Notes:
// - Input is an edge list: each line "u v" (ints). '#' starts a comment. Self-loops skipped.
// - By default the graph is UNDIRECTED.
// - If --directed is passed, degree/density reflect directed arcs,
//   but components, shortest paths, and clustering are computed on the UNDIRECTED projection (u<->v).
//   This keeps small-world style metrics intuitive.
//
// Build (Apple clang + Homebrew libomp):
//   brew install libomp
//   clang -O3 -std=c11 -Xpreprocessor -fopenmp \
//     -I"$(brew --prefix libomp)/include" graph_metrics_omp.c \
//     -L"$(brew --prefix libomp)/lib" -lomp -o graph_metrics_omp -mcpu=native
//
// Usage:
//   ./graph_metrics_omp edge_list.txt [--per-node] [--directed] [--threads N]

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <omp.h>

typedef struct { int *data; int size; int cap; } IntVec;
typedef struct { int a,b; } Edge;
typedef struct { Edge *data; int size; int cap; } EdgeVec;

static void ivec_init(IntVec *v){ v->data=NULL; v->size=0; v->cap=0; }
static void ivec_free(IntVec *v){ free(v->data); v->data=NULL; v->size=v->cap=0; }
static void ivec_push(IntVec *v, int x){
    if(v->size==v->cap){ v->cap = v->cap? v->cap*2 : 4; v->data=(int*)realloc(v->data, v->cap*sizeof(int)); }
    v->data[v->size++] = x;
}
static void evec_init(EdgeVec *v){ v->data=NULL; v->size=0; v->cap=0; }
static void evec_free(EdgeVec *v){ free(v->data); v->data=NULL; v->size=v->cap=0; }
static void evec_push(EdgeVec *v, int a, int b){
    if(v->size==v->cap){ v->cap=v->cap? v->cap*2:4; v->data=(Edge*)realloc(v->data, v->cap*sizeof(Edge)); }
    v->data[v->size].a=a; v->data[v->size].b=b; v->size++;
}

static int cmp_int(const void*pa,const void*pb){
    int a=*(const int*)pa, b=*(const int*)pb;
    return (a<b)?-1:(a>b);
}
static int cmp_edge(const void*pa,const void*pb){
    const Edge *A=(const Edge*)pa, *B=(const Edge*)pb;
    if(A->a!=B->a) return (A->a<B->a)?-1:1;
    if(A->b!=B->b) return (A->b<B->b)?-1:1;
    return 0;
}
static int unique_ints(int *a,int n){
    if(n==0) return 0;
    int j=1;
    for(int i=1;i<n;i++) if(a[i]!=a[i-1]) a[j++]=a[i];
    return j;
}
static int unique_edges(Edge *e,int n){
    if(n==0) return 0;
    int j=1;
    for(int i=1;i<n;i++){
        if(e[i].a!=e[i-1].a || e[i].b!=e[i-1].b) e[j++]=e[i];
    }
    return j;
}

typedef struct {
    int n;                 // nodes
    long long m_dir;       // directed arcs (if directed) or 2*undirected edges (if you want)
    long long m_undir;     // undirected edges count
    int directed;          // input flag
    IntVec *adj;           // directed adjacency (out-neighbors) if directed, else same as uadj
    IntVec *uadj;          // undirected adjacency (neighbors both ways), used for components/APSP/clustering
    int *orig_id;          // compact index -> original id
    int *deg_out;          // out-degree (if directed), equals degree if undirected
    int *deg_und;          // undirected degree
} Graph;

// reading helpers
static int *raw_ids=NULL; static int raw_ids_sz=0, raw_ids_cap=0;
static void push_id(int x){
    if(raw_ids_sz==raw_ids_cap){
        raw_ids_cap = raw_ids_cap? raw_ids_cap*2 : 1024;
        raw_ids = (int*)realloc(raw_ids, raw_ids_cap*sizeof(int));
    }
    raw_ids[raw_ids_sz++]=x;
}
static int is_comment_or_blank(const char*s){
    while(*s && isspace((unsigned char)*s)) s++;
    return (*s=='#' || *s=='\0');
}

// global map of unique sorted ids for compacting
static int *uniq_ids=NULL; static int uniq_n=0;
static int id_to_idx(int id){
    int lo=0, hi=uniq_n-1;
    while(lo<=hi){
        int mid=(lo+hi)>>1;
        int v=uniq_ids[mid];
        if(v==id) return mid;
        if(v<id) lo=mid+1; else hi=mid-1;
    }
    return -1; // should not happen
}

static void free_graph(Graph *g){
    if(!g) return;
    if(g->adj){ for(int i=0;i<g->n;i++) ivec_free(&g->adj[i]); free(g->adj); }
    if(g->uadj){ for(int i=0;i<g->n;i++) ivec_free(&g->uadj[i]); free(g->uadj); }
    free(g->deg_out); free(g->deg_und);
    free(g->orig_id);
    free(raw_ids); raw_ids=NULL; raw_ids_sz=raw_ids_cap=0;
    free(uniq_ids); uniq_ids=NULL; uniq_n=0;
    memset(g,0,sizeof(*g));
}

static int build_graph(const char*fname, int directed, Graph *g){
    memset(g,0,sizeof(*g));
    g->directed = directed?1:0;

    FILE *fp=fopen(fname,"r");
    if(!fp){ fprintf(stderr,"Error: cannot open %s\n", fname); return 0; }

    EdgeVec E; evec_init(&E);
    char buf[1<<16];
    while(fgets(buf,sizeof(buf),fp)){
        if(is_comment_or_blank(buf)) continue;
        int u,v; if(sscanf(buf,"%d %d",&u,&v)!=2) continue;
        if(u==v) continue; // skip self-loops
        if(!g->directed && u>v){ int t=u; u=v; v=t; } // normalize for undirected
        evec_push(&E,u,v);
        push_id(u); push_id(v);
    }
    fclose(fp);
    if(E.size==0){ fprintf(stderr,"No edges read.\n"); evec_free(&E); return 0; }

    // compact ids
    qsort(raw_ids, raw_ids_sz, sizeof(int), cmp_int);
    uniq_n = unique_ints(raw_ids, raw_ids_sz);
    uniq_ids=(int*)malloc(uniq_n*sizeof(int));
    for(int i=0;i<uniq_n;i++) uniq_ids[i]=raw_ids[i];

    g->n = uniq_n;
    g->orig_id=(int*)malloc(g->n*sizeof(int));
    for(int i=0;i<g->n;i++) g->orig_id[i]=uniq_ids[i];

    // map edges to compact indices
    for(int i=0;i<E.size;i++){
        E.data[i].a = id_to_idx(E.data[i].a);
        E.data[i].b = id_to_idx(E.data[i].b);
    }
    // dedupe edges
    if(!g->directed){
        qsort(E.data, E.size, sizeof(Edge), cmp_edge);
        E.size = unique_edges(E.data, E.size);
    }else{
        qsort(E.data, E.size, sizeof(Edge), cmp_edge);
        E.size = unique_edges(E.data, E.size); // dedupe identical arcs
    }

    // allocate adjacencies
    g->adj  = (IntVec*)malloc(g->n*sizeof(IntVec));
    g->uadj = (IntVec*)malloc(g->n*sizeof(IntVec));
    g->deg_out=(int*)calloc(g->n,sizeof(int));
    g->deg_und=(int*)calloc(g->n,sizeof(int));
    for(int i=0;i<g->n;i++){ ivec_init(&g->adj[i]); ivec_init(&g->uadj[i]); }

    // fill adjacencies
    long long m_undir = 0, m_dir = 0;
    if(!g->directed){
        for(int i=0;i<E.size;i++){
            int u=E.data[i].a, v=E.data[i].b;
            ivec_push(&g->adj[u], v);
            ivec_push(&g->adj[v], u);
            ivec_push(&g->uadj[u], v);
            ivec_push(&g->uadj[v], u);
            m_undir++;
        }
        m_dir = 2*m_undir; // conceptual
    }else{
        // directed arcs
        for(int i=0;i<E.size;i++){
            int u=E.data[i].a, v=E.data[i].b;
            ivec_push(&g->adj[u], v);
            // undirected projection (for components/APSP/clustering)
            ivec_push(&g->uadj[u], v);
            ivec_push(&g->uadj[v], u);
            m_dir++;
        }
        // count undirected unique edges from projection (approximate via set later after sort+unique)
    }

    // sort+unique adj lists; compute degrees and finalize m_undir for directed case
    long long undir_sum = 0;
    for(int i=0;i<g->n;i++){
        // directed list
        if(g->adj[i].size>1){
            qsort(g->adj[i].data, g->adj[i].size, sizeof(int), cmp_int);
            g->adj[i].size = unique_ints(g->adj[i].data, g->adj[i].size);
        }
        g->deg_out[i] = g->adj[i].size;
        // undirected list
        if(g->uadj[i].size>1){
            qsort(g->uadj[i].data, g->uadj[i].size, sizeof(int), cmp_int);
            g->uadj[i].size = unique_ints(g->uadj[i].data, g->uadj[i].size);
        }
        g->deg_und[i] = g->uadj[i].size;
        undir_sum += g->deg_und[i];
    }
    if(!g->directed){
        g->m_undir = m_undir; g->m_dir = m_dir;
    }else{
        // each undirected edge counted twice in degree sum
        g->m_undir = undir_sum/2; g->m_dir = m_dir;
    }

    evec_free(&E);
    return 1;
}

// binary search membership in sorted neighbor list
static int is_adjacent(const IntVec *adj, int u, int v){
    const int *a = adj[u].data;
    int lo=0, hi=adj[u].size-1;
    while(lo<=hi){
        int mid=(lo+hi)>>1;
        int w=a[mid];
        if(w==v) return 1;
        if(w<v) lo=mid+1; else hi=mid-1;
    }
    return 0;
}

// BFS over UNDIRECTED adjacency (for components/APSP)
static int bfs_component(const Graph *g, int s, int comp_id, int *comp_of, int *queue){
    int head=0, tail=0, size=0;
    queue[tail++]=s; comp_of[s]=comp_id;
    while(head<tail){
        int u=queue[head++]; size++;
        const IntVec *U = &g->uadj[u];
        for(int j=0;j<U->size;j++){
            int v=U->data[j];
            if(comp_of[v]==-1){
                comp_of[v]=comp_id;
                queue[tail++]=v;
            }
        }
    }
    return size;
}

// Count neighbor-neighbor links for node v using UNDIRECTED adjacency
static long long neighbor_links_undirected(const Graph *g, int v){
    const int *nbr = g->uadj[v].data;
    int k = g->uadj[v].size;
    if(k<2) return 0;
    long long links=0;
    for(int i=0;i<k;i++){
        int a=nbr[i];
        for(int j=i+1;j<k;j++){
            int b=nbr[j];
            if(is_adjacent(g->uadj, a, b)) links++;
        }
    }
    return links;
}

// APSP (undirected projection): parallel over sources
static void apsp_avg_and_diam_parallel(const Graph *g, double *avg_out, int *diam_out){
    int n=g->n;
    long long global_pairs=0;
    long long global_sum=0;
    int global_diam=0;

    #pragma omp parallel
    {
        int *dist = (int*)malloc(n*sizeof(int));
        int *queue= (int*)malloc(n*sizeof(int));

        long long local_pairs=0, local_sum=0;
        int local_diam=0;

        #pragma omp for schedule(dynamic,1) nowait
        for(int s=0; s<n; s++){
            for(int i=0;i<n;i++) dist[i]=-1;
            int head=0, tail=0;
            dist[s]=0; queue[tail++]=s;
            int far=0;

            while(head<tail){
                int u=queue[head++];
                const IntVec *U = &g->uadj[u];
                for(int j=0;j<U->size;j++){
                    int v=U->data[j];
                    if(dist[v]==-1){
                        dist[v]=dist[u]+1;
                        if(dist[v]>far) far=dist[v];
                        queue[tail++]=v;
                    }
                }
            }
            if(far>local_diam) local_diam=far;
            for(int t=s+1;t<n;t++){
                if(dist[t]>=0){ local_sum += dist[t]; local_pairs++; }
            }
        }

        #pragma omp atomic
        global_sum += local_sum;
        #pragma omp atomic
        global_pairs += local_pairs;
        #pragma omp critical
        { if(local_diam>global_diam) global_diam=local_diam; }

        free(queue); free(dist);
    }

    *diam_out = global_diam;
    *avg_out = (global_pairs>0)? (double)global_sum/(double)global_pairs : 0.0;
}

int main(int argc, char **argv){
    if(argc<2){
        fprintf(stderr,"Usage: %s edge_list.txt [--per-node] [--directed] [--threads N]\n", argv[0]);
        return 1;
    }
    const char *fname = argv[1];
    int want_per_node = 0;
    int directed = 0;
    int threads = 0;

    for(int i=2;i<argc;i++){
        if(strcmp(argv[i],"--per-node")==0) want_per_node=1;
        else if(strcmp(argv[i],"--directed")==0) directed=1;
        else if(strcmp(argv[i],"--threads")==0 && i+1<argc){
            threads = atoi(argv[++i]);
            if(threads<1) threads=1;
        }
    }

    if(threads>0){
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
    }

    Graph G;
    if(!build_graph(fname, directed, &G)) return 1;

    // Basic size metrics (undirected degree & density are primary)
    long long deg_sum_und=0, deg_sum_out=0;
    int min_deg = (G.n>0)? G.deg_und[0]:0;
    int max_deg = 0;
    for(int i=0;i<G.n;i++){
        deg_sum_und += G.deg_und[i];
        deg_sum_out += G.deg_out[i];
        if(G.deg_und[i]<min_deg) min_deg=G.deg_und[i];
        if(G.deg_und[i]>max_deg) max_deg=G.deg_und[i];
    }
    double avg_deg_und = (G.n>0)? (double)deg_sum_und/(double)G.n : 0.0;
    double avg_deg_out = (G.n>0)? (double)deg_sum_out/(double)G.n : 0.0;
    double density = 0.0;
    if(!G.directed){
        density = (G.n>1)? (2.0*(double)G.m_undir)/((double)G.n*(G.n-1)) : 0.0;
    }else{
        density = (G.n>1)? ((double)G.m_dir)/((double)G.n*(G.n-1)) : 0.0;
    }

    // Connected components (weak; on undirected projection)
    int *comp_of=(int*)malloc(G.n*sizeof(int));
    for(int i=0;i<G.n;i++) comp_of[i]=-1;
    int *queue=(int*)malloc(G.n*sizeof(int));
    int comp_cnt=0;
    // collect sizes
    IntVec comp_sizes; ivec_init(&comp_sizes);

    for(int i=0;i<G.n;i++){
        if(comp_of[i]==-1){
            int sz=bfs_component(&G, i, comp_cnt, comp_of, queue);
            ivec_push(&comp_sizes, sz);
            comp_cnt++;
        }
    }
    // sort sizes desc
    if(comp_sizes.size>1){
        qsort(comp_sizes.data, comp_sizes.size, sizeof(int), cmp_int);
        for(int i=0,j=comp_sizes.size-1;i<j;i++,j--){
            int t=comp_sizes.data[i]; comp_sizes.data[i]=comp_sizes.data[j]; comp_sizes.data[j]=t;
        }
    }

    // Clustering & triangles (UNDIRECTED) — parallelized per node
    double avg_local_C = 0.0;
    long long total_triangles = 0;
    {
        double sum_C=0.0;
        long long sum_tri_touch=0;

        #pragma omp parallel for schedule(static) reduction(+:sum_C,sum_tri_touch)
        for(int v=0; v<G.n; v++){
            long long tri_v = neighbor_links_undirected(&G, v);
            int k = G.deg_und[v];
            double Cv = 0.0;
            if(k>=2) Cv = (2.0*(double)tri_v)/((double)k*(k-1));
            sum_C += Cv;
            sum_tri_touch += tri_v;
        }
        avg_local_C = (G.n>0)? sum_C/(double)G.n : 0.0;
        total_triangles = sum_tri_touch/3;
    }

    // APSP average and diameter on UNDIRECTED projection — parallel
    double avg_shortest=0.0; int diameter=0;
    apsp_avg_and_diam_parallel(&G, &avg_shortest, &diameter);

    // -------- Output --------
    printf("=== Graph Summary ===\n");
    printf("Directed input: %s\n", G.directed? "Yes":"No");
    printf("Nodes (n): %d\n", G.n);
    if(!G.directed){
        printf("Edges (m, undirected): %lld\n", G.m_undir);
    }else{
        printf("Arcs (m, directed): %lld\n", G.m_dir);
        printf("Undirected projection edges: %lld\n", G.m_undir);
    }
    printf("Average degree (undirected view): %.6f\n", avg_deg_und);
    if(G.directed) printf("Average out-degree (directed): %.6f\n", avg_deg_out);
    printf("Min degree (undirected): %d, Max degree (undirected): %d\n", min_deg, max_deg);
    printf("Density: %.8f\n", density);

    printf("\n=== Connectivity (weak, undirected projection) ===\n");
    printf("Connected components: %d\n", comp_cnt);
    printf("Component sizes (desc, top 10):");
    int lim = comp_sizes.size<10? comp_sizes.size:10;
    for(int i=0;i<lim;i++) printf(" %d", comp_sizes.data[i]);
    if(comp_sizes.size>10) printf(" ...");
    printf("\n");

    printf("\n=== Shortest Paths (undirected projection) ===\n");
    printf("Average shortest path length (reachable pairs): %.6f\n", avg_shortest);
    printf("Diameter (longest finite shortest path): %d\n", diameter);

    printf("\n=== Clustering & Triangles (undirected) ===\n");
    printf("Average local clustering coefficient: %.6f\n", avg_local_C);
    printf("Total triangles: %lld\n", total_triangles);

    if(want_per_node){
        printf("\n=== Per-node (index original_id) ===\n");
        for(int v=0; v<G.n; v++){
            long long tri_v = neighbor_links_undirected(&G, v);
            int k = G.deg_und[v];
            double Cv = (k>=2)? (2.0*(double)tri_v)/((double)k*(k-1)) : 0.0;
            printf("v=%d (%d): deg=%d, C=%.6f\n", v, G.orig_id[v], G.deg_und[v], Cv);
        }
    }

    // cleanup
    ivec_free(&comp_sizes);
    free(comp_of); free(queue);
    free_graph(&G);
    return 0;
}
