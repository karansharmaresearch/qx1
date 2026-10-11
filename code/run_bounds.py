import json, os, platform, sys, time
from datetime import datetime, timezone

import gmpy2, mpmath

sys.path.insert(0, os.path.dirname(__file__))
import bounds as b

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "data", "bounds")
QS = [3, 5, 7, 15, 31, 63, 127, 255, 511, 1023]


def run(q, Cmax, log):
    t0 = time.time()
    J, Bup, gap = b.certified_bounds(q, Cmax)
    recs = b.records(Bup)
    Llo, Lhi, l2lo = b.log_enclosures(q)
    cf = b.continued_fraction_log2(q, 40)
    conv = b.convergents(cf)
    nchk = b.check_convergents_exact(q, conv, smax=Cmax)
    k = q.bit_length() if (q + 1) & q == 0 else None
    sub = sorted({pow(2, i, q) for i in range(q.bit_length() + 1)}) if q > 1 else []
    res = {
        "q": q, "k": k, "Cmax": Cmax,
        "subgroup_2_mod_q": sub,
        "log2q_enclosure": {"P": b.P, "lo_num": str(Llo), "hi_num": str(Lhi),
                            "note": "lo_num/2^P <= log2 q <= hi_num/2^P"},
        "ln2_lower": {"P": b.P, "num": str(l2lo)},
        "cf_log2q": cf,
        "convergents": [[p, s] for p, s in conv],
        "convergents_exact_checked": nchk,
        "records": [{"c": c, "j": J[c], "Bup": B} for c, B in recs],
        "theorem_table": b.theorem_table(Bup, recs),
        "methodA_vs_B_max_gap_at_records": gap,
        "seconds": round(time.time() - t0, 2),
    }
    with open(os.path.join(OUT, f"q{q}.json"), "w") as f:
        json.dump(res, f, indent=1)
    log(f"q={q:5d}  Cmax={Cmax}  records={len(recs)}  record gap A/B={gap}  "
        f"CF={cf[:12]}...  conv. exact-checked={nchk}  {res['seconds']}s")
    for row in res["theorem_table"][-6:]:
        log(f"    no cycle with c <= {row['C']:>8,}  needs all minima m < {row['B']:>16,} "
            f"(record at c={row['record_c']})")
    return res


def main():
    Cmax = int(sys.argv[1]) if len(sys.argv) > 1 else 500_000
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")

    meta = {"started": datetime.now(timezone.utc).isoformat(),
            "python": sys.version.split()[0], "platform": platform.platform(),
            "gmpy2": gmpy2.version(), "mpmath": mpmath.__version__,
            "P": b.P, "IV_PREC": b.IV_PREC, "Cmax": Cmax}
    log(f"=== B1 run {meta['started']}  {json.dumps(meta)}")
    summary = {"meta": meta, "q": {}}
    for q in QS:
        r = run(q, Cmax, log)
        summary["q"][str(q)] = {"theorem_table": r["theorem_table"],
                                "cf_log2q": r["cf_log2q"][:20]}
    meta["finished"] = datetime.now(timezone.utc).isoformat()
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    log(f"=== finished {meta['finished']}")


if __name__ == "__main__":
    main()
