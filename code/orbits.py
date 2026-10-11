import json, math, os, sys, time
from datetime import datetime, timezone

import gmpy2
from gmpy2 import mpz
from scipy.stats import chisquare

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "data", "orbits")
QS = [7, 15, 31, 63, 127, 255, 511, 1023]
BITS = 2000
VBINS = 12


def step(q, n):
    t = q * n + 1
    v = gmpy2.bit_scan1(t)
    return t >> v, v


def d2_orbit_of_3(q, S):
    n, vs, trace = mpz(3), [0] * (VBINS + 1), [(0, 2)]
    reached_one = None
    for s in range(1, S + 1):
        n, v = step(q, n)
        vs[min(v, VBINS)] += 1
        if n == 1:
            reached_one = s
            break
        if s % 500 == 0:
            trace.append((s, int(n.bit_length())))
    bits = int(n.bit_length())
    steps = reached_one or S
    return {"steps": steps, "reached_one_at": reached_one, "final_bits": bits,
            "bits_per_step": (bits - 2) / steps, "theory_bits_per_step": math.log2(q / 4),
            "trace": trace, "v_counts": vs}


def d1_d3_small_orbits(q, nstart):
    vs = [0] * (VBINS + 1)
    owner, parent = {}, {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    reach1, cert_len = set(), {}
    starts = list(range(1, 2 * nstart, 2))
    for s0 in starts:
        parent[s0] = s0
        n, steps = mpz(s0), 0
        while True:
            if n == 1:
                reach1.add(s0); cert_len[s0] = steps
                break
            key = int(n)
            if key in owner and owner[key] != s0:
                ra, rb = find(owner[key]), find(s0)
                parent[rb] = ra
                break
            owner[key] = s0
            if n.bit_length() > BITS:
                break
            n, v = step(q, n)
            steps += 1
            if steps <= 1000:
                vs[min(v, VBINS)] += 1
    roots1 = {find(s) for s in reach1}
    r1 = sorted(s for s in starts if find(s) in roots1)
    rest = [s for s in starts if find(s) not in roots1]
    comps = {}
    for s in rest:
        comps.setdefault(find(s), []).append(s)
    reps = sorted(min(c) for c in comps.values())
    N = 2 * nstart - 1
    all_n_reaching = sum(1 for o in r1 for e in range(64) if (o << e) <= N)
    return vs, {"N": N, "odd_starts": len(starts), "odd_reach_1": r1,
                "all_n_reach_1_count": all_n_reaching,
                "independent_escaping_orbits": len(comps), "representatives": reps,
                "largest_families": sorted(((len(c), min(c)) for c in comps.values()), reverse=True)[:8],
                "cutoff_bits": BITS}


def chi2(vcounts):
    n = sum(vcounts[1:])
    exp = [n * 2.0 ** -v for v in range(1, VBINS)] + [n * 2.0 ** -(VBINS - 1)]
    obs = vcounts[1:VBINS] + [vcounts[VBINS]]
    st, p = chisquare(obs, exp)
    mean = sum(min(v, VBINS) * c for v, c in enumerate(vcounts)) / n
    return {"n": n, "observed": obs, "expected": exp, "chi2": st, "dof": len(obs) - 1,
            "p_value": p, "mean_v_capped": mean}


def main():
    S = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
    nstart = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    os.makedirs(OUT, exist_ok=True)
    log = open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8")
    log.write(f"=== {datetime.now(timezone.utc).isoformat()} S={S} NSTART={nstart} gmpy2={gmpy2.version()}\n")
    for q in QS:
        t0 = time.time()
        o3 = d2_orbit_of_3(q, S)
        vs_small, merge = d1_d3_small_orbits(q, nstart)
        res = {"q": q, "k": q.bit_length(), "S": S, "orbit_of_3": o3, "merge": merge,
               "v_chi2_orbit3": chi2(o3["v_counts"]), "v_chi2_small": chi2(vs_small),
               "seconds": round(time.time() - t0, 2)}
        json.dump(res, open(os.path.join(OUT, f"q{q}.json"), "w"), indent=1)
        line = (f"q={q}: orbit of 3 -> {o3['final_bits']} bits after {o3['steps']} odd steps "
                f"({o3['bits_per_step']:.4f}/step vs log2(q/4)={o3['theory_bits_per_step']:.4f}); "
                f"chi2(v) orbit3 p={res['v_chi2_orbit3']['p_value']:.3g}, small p={res['v_chi2_small']['p_value']:.3g}; "
                f"odd n<={merge['N']}: {len(merge['odd_reach_1'])} reach 1 (all n: {merge['all_n_reach_1_count']}), "
                f"{merge['independent_escaping_orbits']} independent escaping orbits; {res['seconds']}s")
        print(line, flush=True); log.write(line + "\n")


if __name__ == "__main__":
    main()
