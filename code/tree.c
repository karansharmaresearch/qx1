#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned __int128 u128;
#define MAXD 38

static uint64_t Q;
static int K, D;
static u128 LMAX, LIM;
static u128 P10[MAXD + 1];
static int *DLOG;

static void die(const char *s) { fprintf(stderr, "error: %s\n", s); exit(1); }

static inline int dec(u128 x) { int a = 0; while (P10[a] < x) a++; return a; }

typedef struct { u128 m; int pmd; } node_t;
typedef struct { uint64_t h[MAXD + 1][MAXD + 1]; uint64_t nodes; } hist_t;

#define FOR_CHILDREN(m, pmd, BODY)                                              \
    do {                                                                       \
        int a_ = DLOG[(uint64_t)((m) % Q)];                                    \
        if (a_ >= 0) {                                                         \
            int v0_ = (K - a_ % K) % K; if (v0_ == 0) v0_ = K;                 \
            if ((m) <= (LIM >> v0_)) {                                         \
                u128 cur_ = (m) << v0_;                                        \
                for (;;) {                                                     \
                    u128 n_ = (cur_ - 1) / Q;                                  \
                    if (n_ != 1) {                                             \
                        int dn_ = dec(n_);                                     \
                        int cp_ = dn_ > (pmd) ? dn_ : (pmd);                   \
                        BODY                                                   \
                    }                                                          \
                    if (cur_ > (LIM >> K)) break;                              \
                    cur_ <<= K;                                                \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)

static void dfs(node_t root, hist_t *H) {
    size_t cap = 1 << 16, sp = 0;
    node_t *st = malloc(cap * sizeof *st);
    st[sp++] = root;
    while (sp) {
        node_t x = st[--sp];
        FOR_CHILDREN(x.m, x.pmd, {
            H->h[dn_][cp_]++; H->nodes++;
            if (sp == cap) { cap *= 2; st = realloc(st, cap * sizeof *st); }
            st[sp].m = n_; st[sp].pmd = cp_; sp++;
        });
    }
    free(st);
}

static void print_u128(FILE *f, u128 x) {
    char b[48]; int i = 47; b[i] = 0;
    do { b[--i] = '0' + (int)(x % 10); x /= 10; } while (x);
    fputs(b + i, f);
}

int main(int argc, char **argv) {
    if (argc != 4) die("usage: tree q D OUT.json");
    Q = strtoull(argv[1], 0, 10); D = atoi(argv[2]);
    if (Q < 3 || (Q & (Q + 1))) die("q must be 2^k - 1 >= 3");
    K = 0; while ((1ULL << K) <= Q) K++;
    if (D < 1 || D > 37) die("D out of range");
    P10[0] = 1; for (int i = 1; i <= MAXD; i++) P10[i] = P10[i - 1] * 10;
    LMAX = P10[D];
    if (LMAX > (((u128)1 << 126) / Q)) die("q*10^D too large for 128-bit arithmetic");
    LIM = (u128)Q * LMAX + 1;
    DLOG = malloc(Q * sizeof *DLOG);
    for (uint64_t r = 0; r < Q; r++) DLOG[r] = -1;
    { uint64_t x = 1; for (int a = 0; a < K; a++) { DLOG[x] = a; x = 2 * x % Q; } }

    double t0 = omp_get_wtime();
    hist_t *tot = calloc(1, sizeof *tot);
    tot->h[0][0] = 1; tot->nodes = 1;

    size_t fcap = 1024, fn = 0;
    node_t *fr = malloc(fcap * sizeof *fr);
    fr[fn++] = (node_t){1, 0};
    int levels = 0;
    while (fn > 0 && fn < 200000 && levels < 64) {
        size_t ncap = 1024, nn = 0;
        node_t *nx = malloc(ncap * sizeof *nx);
        for (size_t i = 0; i < fn; i++) {
            FOR_CHILDREN(fr[i].m, fr[i].pmd, {
                tot->h[dn_][cp_]++; tot->nodes++;
                if (nn == ncap) { ncap *= 2; nx = realloc(nx, ncap * sizeof *nx); }
                nx[nn].m = n_; nx[nn].pmd = cp_; nn++;
            });
        }
        free(fr); fr = nx; fn = nn; levels++;
    }
    #pragma omp parallel
    {
        hist_t *loc = calloc(1, sizeof *loc);
        #pragma omp for schedule(dynamic, 64)
        for (int64_t i = 0; i < (int64_t)fn; i++) dfs(fr[i], loc);
        #pragma omp critical
        {
            for (int a = 0; a <= MAXD; a++)
                for (int b = 0; b <= MAXD; b++) tot->h[a][b] += loc->h[a][b];
            tot->nodes += loc->nodes;
        }
        free(loc);
    }
    double sec = omp_get_wtime() - t0;

    FILE *f = fopen(argv[3], "w");
    if (!f) die("cannot write output");
    fprintf(f, "{\n \"q\": %llu, \"k\": %d, \"D\": %d, \"Lmax\": \"", (unsigned long long)Q, K, D);
    print_u128(f, LMAX);
    fprintf(f, "\",\n \"nodes\": %llu, \"bfs_levels\": %d, \"frontier\": %llu,\n",
            (unsigned long long)tot->nodes, levels, (unsigned long long)fn);
    fprintf(f, " \"seconds\": %.3f, \"threads\": %d, \"compiler\": \"%s\",\n", sec,
            omp_get_max_threads(), __VERSION__);
    fprintf(f, " \"hist_note\": \"H[a][b] = #{n: dec(n)=a, dec(orbit max)=b}, dec(x)=min{a: x<=10^a}\",\n");
    fprintf(f, " \"H\": [\n");
    for (int a = 0; a <= D; a++) {
        fprintf(f, "  [");
        for (int b = 0; b <= D; b++) fprintf(f, "%llu%s", (unsigned long long)tot->h[a][b], b < D ? ", " : "");
        fprintf(f, "]%s\n", a < D ? "," : "");
    }
    fprintf(f, " ]\n}\n");
    fclose(f);
    printf("q=%llu D=%d nodes=%llu %.2fs\n", (unsigned long long)Q, D, (unsigned long long)tot->nodes, sec);
    return 0;
}
