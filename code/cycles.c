#include <gmp.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned __int128 u128;

static uint64_t Q, CMAXALL, BOUND;
static int64_t NST;
static uint64_t *ST_C, *ST_B;
static unsigned char *ADM;

static void die(const char *msg) { fprintf(stderr, "error: %s\n", msg); exit(1); }

static void read_staircase(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) die("cannot open staircase file");
    char line[1 << 16];
    int64_t cap = 0;
    NST = 0;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#') {
            char key[64];
            if (sscanf(line, "# %63s", key) != 1) continue;
            if (!strcmp(key, "q")) sscanf(line, "# q %llu", (unsigned long long *)&Q);
            else if (!strcmp(key, "C")) sscanf(line, "# C %llu", (unsigned long long *)&CMAXALL);
            else if (!strcmp(key, "B")) sscanf(line, "# B %llu", (unsigned long long *)&BOUND);
            else if (!strcmp(key, "subgroup")) {
                if (!Q) die("q must precede subgroup");
                ADM = calloc(Q, 1);
                char *p = line + strlen("# subgroup");
                char *end;
                for (;;) {
                    unsigned long long r = strtoull(p, &end, 10);
                    if (end == p) break;
                    if (r >= Q) die("subgroup residue >= q");
                    ADM[r] = 1;
                    p = end;
                }
            }
            continue;
        }
        unsigned long long c, B;
        if (sscanf(line, "%llu %llu", &c, &B) != 2) continue;
        if (NST == cap) {
            cap = cap ? 2 * cap : 1024;
            ST_C = realloc(ST_C, cap * sizeof *ST_C);
            ST_B = realloc(ST_B, cap * sizeof *ST_B);
        }
        ST_C[NST] = c; ST_B[NST] = B; NST++;
    }
    fclose(f);
    if (!Q || !CMAXALL || !BOUND || !ADM || NST == 0) die("incomplete staircase file");
    for (int64_t i = 1; i < NST; i++)
        if (!(ST_C[i] > ST_C[i - 1] && ST_B[i] < ST_B[i - 1])) die("staircase not monotone");
    if (ST_B[0] != BOUND) die("B != first staircase bound");
    if (ST_C[NST - 1] != CMAXALL) die("last staircase c != C");
}

static uint64_t cmax_of(uint64_t m) {
    int64_t lo = 0, hi = NST;
    while (lo < hi) {
        int64_t mid = (lo + hi) / 2;
        if (ST_B[mid] > m) lo = mid + 1; else hi = mid;
    }
    return lo ? ST_C[lo - 1] : 0;
}

enum { OUT_DROP = 0, OUT_SURV = 1, OUT_FOUND = 2 };

static inline int ctz128(u128 x) {
    uint64_t lo = (uint64_t)x;
    return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((uint64_t)(x >> 64));
}

static void mpz_set_u128(mpz_t z, u128 x) {
    uint64_t w[2] = {(uint64_t)x, (uint64_t)(x >> 64)};
    mpz_import(z, 2, -1, sizeof(uint64_t), 0, 0, w);
}

static u128 mpz_get_u128(const mpz_t z) {
    uint64_t w[2] = {0, 0};
    mpz_export(w, NULL, -1, sizeof(uint64_t), 0, 0, z);
    return ((u128)w[1] << 64) | w[0];
}

static int BLOCK = 1;

static uint64_t block_step(mpz_t n, uint64_t maxsteps, mpz_t az, mpz_t rz) {
    u128 y = (u128)mpz_getlimbn(n, 0);
    u128 r = 0;
    uint64_t A = 1, i = 0;
    int S = 0, p = 64;
    const uint64_t ALIM = ((uint64_t)1 << 56) / Q;
    while (i < maxsteps && A <= ALIM && p >= 8) {
        u128 t = ((u128)Q * y + 1) & ((((u128)1) << p) - 1);
        if (t == 0) break;
        int v = ctz128(t);
        r = (u128)Q * r + ((u128)1 << S);
        A *= Q; S += v; p -= v;
        y = t >> v;
        i++;
    }
    if (i == 0) return 0;
    mpz_import(az, 1, -1, sizeof A, 0, 0, &A);
    mpz_set_u128(rz, r);
    mpz_mul(n, n, az);
    mpz_add(n, n, rz);
    if (mpz_scan1(n, 0) != (mp_bitcnt_t)S) die("internal: block step not exact");
    mpz_tdiv_q_2exp(n, n, S);
    return i;
}

static int run_m(uint64_t m, uint64_t K, mpz_t n, mpz_t mz, mpz_t az, mpz_t rz,
                 uint64_t *step) {
    const u128 LIM = ((~(u128)0) - 1) / Q;
    const size_t mbits = 64 - __builtin_clzll(m);
    u128 x = m;
    int big = 0;
    uint64_t s = 0;
    while (s < K) {
        if (!big) {
            if (x <= LIM) {
                u128 t = (u128)Q * x + 1;
                x = t >> ctz128(t);
                s++;
                if (x < m) { *step = s; return OUT_DROP; }
                if (x == m) { *step = s; return OUT_FOUND; }
                continue;
            }
            mpz_set_u128(n, x);
            big = 1;
        }
        if (BLOCK && K - s >= 2 && mpz_sizeinbase(n, 2) > mbits + 66) {
            uint64_t i = block_step(n, K - s, az, rz);
            if (i) {
                s += i;
                if (mpz_cmp(n, mz) <= 0) die("internal: block reached <= m");
                continue;
            }
        }
        mpz_mul_ui(n, n, (unsigned long)Q);
        mpz_add_ui(n, n, 1);
        mpz_tdiv_q_2exp(n, n, mpz_scan1(n, 0));
        s++;
        int cmp = mpz_cmp(n, mz);
        if (cmp < 0) { *step = s; return OUT_DROP; }
        if (cmp == 0) { *step = s; return OUT_FOUND; }
        if (mpz_sizeinbase(n, 2) <= 100) {
            x = mpz_get_u128(n);
            big = 0;
        }
    }
    *step = K;
    return OUT_SURV;
}

static void mpz_set_u64(mpz_t z, uint64_t v) { mpz_import(z, 1, -1, sizeof v, 0, 0, &v); }

typedef struct { uint64_t cand, drop, surv, found, steps, surv_steps; } stats_t;

static void add(stats_t *a, const stats_t *b) {
    a->cand += b->cand; a->drop += b->drop; a->surv += b->surv;
    a->found += b->found; a->steps += b->steps; a->surv_steps += b->surv_steps;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: cycles STAIRCASE OUTDIR [--mlo M] [--mhi M] [--chunk N] [--dump FILE] [--noblock]\n");
        return 2;
    }
    const char *stfile = argv[1], *outdir = argv[2], *dumpfile = NULL;
    uint64_t mlo = 3, mhi = 0, chunk = (uint64_t)1 << 24;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--mlo") && i + 1 < argc) mlo = strtoull(argv[++i], 0, 10);
        else if (!strcmp(argv[i], "--mhi") && i + 1 < argc) mhi = strtoull(argv[++i], 0, 10);
        else if (!strcmp(argv[i], "--chunk") && i + 1 < argc) chunk = strtoull(argv[++i], 0, 10);
        else if (!strcmp(argv[i], "--dump") && i + 1 < argc) dumpfile = argv[++i];
        else if (!strcmp(argv[i], "--noblock")) BLOCK = 0;
        else die("bad argument");
    }
    if (GMP_LIMB_BITS != 64) die("block steps need 64-bit GMP limbs");
    read_staircase(stfile);
    if (mhi == 0 || mhi > BOUND) mhi = BOUND;
    if (mlo < 3) mlo = 3;
    if (chunk < 2) die("chunk too small");
    if (Q >= ((uint64_t)1 << 27)) die("q too large for the 2^100 fast-path threshold");

    char ckpath[4096], foundpath[4096], respath[4096];
    snprintf(ckpath, sizeof ckpath, "%s/checkpoint_q%llu_C%llu_%llu_%llu.txt", outdir,
             (unsigned long long)Q, (unsigned long long)CMAXALL,
             (unsigned long long)mlo, (unsigned long long)mhi);
    snprintf(foundpath, sizeof foundpath, "%s/found_q%llu_C%llu.txt", outdir,
             (unsigned long long)Q, (unsigned long long)CMAXALL);
    snprintf(respath, sizeof respath, "%s/result_q%llu_C%llu_%llu_%llu.json", outdir,
             (unsigned long long)Q, (unsigned long long)CMAXALL,
             (unsigned long long)mlo, (unsigned long long)mhi);

    stats_t tot = {0};
    uint64_t resume_from = mlo;
    double prev_sec = 0;
    if (!dumpfile) {
        FILE *ck = fopen(ckpath, "r");
        if (ck) {
            char line[512];
            while (fgets(line, sizeof line, ck)) {
                unsigned long long lo, hi, a, b, c, d, e, g; double sec;
                if (line[0] == '#') continue;
                if (sscanf(line, "%llu %llu %llu %llu %llu %llu %llu %llu %lf",
                           &lo, &hi, &a, &b, &c, &d, &e, &g, &sec) == 9) {
                    if (lo != resume_from) die("checkpoint chunks not contiguous");
                    stats_t s = {a, b, c, d, e, g};
                    add(&tot, &s);
                    resume_from = hi; prev_sec += sec;
                }
            }
            fclose(ck);
            fprintf(stderr, "resuming at m = %llu\n", (unsigned long long)resume_from);
        } else {
            ck = fopen(ckpath, "w");
            if (!ck) die("cannot create checkpoint");
            fprintf(ck, "# q %llu C %llu B %llu mlo %llu mhi %llu chunk %llu staircase %s\n",
                    (unsigned long long)Q, (unsigned long long)CMAXALL, (unsigned long long)BOUND,
                    (unsigned long long)mlo, (unsigned long long)mhi,
                    (unsigned long long)chunk, stfile);
            fclose(ck);
        }
    }

    FILE *dump = dumpfile ? fopen(dumpfile, "w") : NULL;
    if (dumpfile && !dump) die("cannot open dump file");
    double t0 = omp_get_wtime();

    for (uint64_t lo = resume_from; lo < mhi; lo += chunk) {
        uint64_t hi = lo + chunk < mhi ? lo + chunk : mhi;
        stats_t cs = {0};
        double tc = omp_get_wtime();
        if (dump) {
            mpz_t n, mz, az, rz; mpz_init(n); mpz_init(mz); mpz_init(az); mpz_init(rz);
            for (uint64_t m = lo | 1; m < hi; m += 2) {
                if (!ADM[m % Q]) continue;
                uint64_t K = cmax_of(m), st;
                if (!K) continue;
                mpz_set_u64(mz, m);
                int o = run_m(m, K, n, mz, az, rz, &st);
                fprintf(dump, "%llu %c %llu\n", (unsigned long long)m, "DSF"[o], (unsigned long long)st);
                cs.cand++; cs.steps += st;
                if (o == OUT_DROP) cs.drop++;
                else if (o == OUT_SURV) { cs.surv++; cs.surv_steps += st; }
                else cs.found++;
            }
            mpz_clear(n); mpz_clear(mz); mpz_clear(az); mpz_clear(rz);
        } else {
            #pragma omp parallel
            {
                mpz_t n, mz, az, rz; mpz_init(n); mpz_init(mz); mpz_init(az); mpz_init(rz);
                stats_t ls = {0};
                int64_t first = (int64_t)(lo | 1), last = (int64_t)hi;
                #pragma omp for schedule(dynamic, 4096) nowait
                for (int64_t mm = first; mm < last; mm += 2) {
                    uint64_t m = (uint64_t)mm;
                    if (!ADM[m % Q]) continue;
                    uint64_t K = cmax_of(m), st;
                    if (!K) continue;
                    mpz_set_u64(mz, m);
                    int o = run_m(m, K, n, mz, az, rz, &st);
                    ls.cand++; ls.steps += st;
                    if (o == OUT_DROP) ls.drop++;
                    else if (o == OUT_SURV) { ls.surv++; ls.surv_steps += st; }
                    else {
                        ls.found++;
                        #pragma omp critical
                        {
                            FILE *ff = fopen(foundpath, "a");
                            if (ff) { fprintf(ff, "%llu %llu\n", (unsigned long long)m, (unsigned long long)st); fclose(ff); }
                            fprintf(stderr, "FOUND cycle: min %llu, %llu odd steps\n",
                                    (unsigned long long)m, (unsigned long long)st);
                        }
                    }
                }
                #pragma omp critical
                add(&cs, &ls);
                mpz_clear(n); mpz_clear(mz); mpz_clear(az); mpz_clear(rz);
            }
        }
        double sec = omp_get_wtime() - tc;
        add(&tot, &cs);
        if (!dump) {
            FILE *ck = fopen(ckpath, "a");
            if (!ck) die("cannot append checkpoint");
            fprintf(ck, "%llu %llu %llu %llu %llu %llu %llu %llu %.3f\n",
                    (unsigned long long)lo, (unsigned long long)hi,
                    (unsigned long long)cs.cand, (unsigned long long)cs.drop,
                    (unsigned long long)cs.surv, (unsigned long long)cs.found,
                    (unsigned long long)cs.steps, (unsigned long long)cs.surv_steps, sec);
            fclose(ck);
        }
        fprintf(stderr, "m in [%llu,%llu): cand %llu drop %llu surv %llu found %llu  %.1fs\n",
                (unsigned long long)lo, (unsigned long long)hi,
                (unsigned long long)cs.cand, (unsigned long long)cs.drop,
                (unsigned long long)cs.surv, (unsigned long long)cs.found, sec);
    }
    double wall = omp_get_wtime() - t0 + prev_sec;
    if (dump) fclose(dump);

    FILE *r = fopen(respath, "w");
    if (!r) die("cannot write result");
    fprintf(r, "{\n \"q\": %llu, \"C\": %llu, \"B\": %llu, \"mlo\": %llu, \"mhi\": %llu,\n",
            (unsigned long long)Q, (unsigned long long)CMAXALL, (unsigned long long)BOUND,
            (unsigned long long)mlo, (unsigned long long)mhi);
    fprintf(r, " \"candidates\": %llu, \"drop\": %llu, \"survive\": %llu, \"found\": %llu,\n",
            (unsigned long long)tot.cand, (unsigned long long)tot.drop,
            (unsigned long long)tot.surv, (unsigned long long)tot.found);
    fprintf(r, " \"total_steps\": %llu, \"survivor_steps\": %llu,\n",
            (unsigned long long)tot.steps, (unsigned long long)tot.surv_steps);
    fprintf(r, " \"seconds\": %.3f, \"threads\": %d, \"gmp\": \"%s\", \"compiler\": \"%s\",\n",
            wall, dump ? 1 : omp_get_max_threads(), gmp_version, __VERSION__);
    fprintf(r, " \"staircase\": \"%s\", \"dump\": %s, \"block\": %s\n}\n", stfile,
            dump ? "true" : "false", BLOCK ? "true" : "false");
    fclose(r);
    printf("q=%llu C=%llu B=%llu m in [%llu,%llu): candidates %llu, drop %llu, survive %llu, FOUND %llu, %.1fs\n",
           (unsigned long long)Q, (unsigned long long)CMAXALL, (unsigned long long)BOUND,
           (unsigned long long)mlo, (unsigned long long)mhi,
           (unsigned long long)tot.cand, (unsigned long long)tot.drop,
           (unsigned long long)tot.surv, (unsigned long long)tot.found, wall);
    return tot.found ? 3 : 0;
}
