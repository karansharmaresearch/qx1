import json, os, subprocess, sys, time
from datetime import datetime, timezone

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "code"))
import bounds as b

OUT = os.path.join(ROOT, "data", "checks")
EXE = os.path.join(ROOT, "code", "cycles.exe")
MAKE_ST = os.path.join(ROOT, "code", "make_staircase.py")


def staircase_path(q, C):
    p = os.path.join(ROOT, "data", "bounds", f"staircase_q{q}_C{C}.txt")
    if not os.path.exists(p):
        subprocess.run([sys.executable, MAKE_ST, str(q), str(C)], check=True,
                       stdout=subprocess.DEVNULL)
    return p


def run_c(q, C, mhi=None, dump=None, noblock=False):
    args = [EXE, staircase_path(q, C), OUT] + (["--noblock"] if noblock else [])
    if mhi:
        args += ["--mhi", str(mhi)]
    if dump:
        args += ["--dump", dump]
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode not in (0, 3):
        raise RuntimeError(r.stderr)
    return r.stdout.strip()


def python_reference(q, C, mhi=None):
    _, BA = b.bounds_exact(q, C)
    B = max(BA[1:])
    mhi = min(mhi or B, B)
    adm = set()
    x = 1
    for _ in range(q):
        adm.add(x)
        x = 2 * x % q
    out = {}
    m_values = [m for m in range(3, mhi, 2) if m % q in adm]
    import bisect
    pairs = sorted((BA[c], c) for c in range(1, C + 1))
    suf = [0] * (len(pairs) + 1)
    for i in range(len(pairs) - 1, -1, -1):
        suf[i] = max(suf[i + 1], pairs[i][1])
    keys = [p[0] for p in pairs]
    for m in m_values:
        K = suf[bisect.bisect_right(keys, m)]
        if K == 0:
            continue
        n, res = m, ("S", K)
        for s in range(1, K + 1):
            t = q * n + 1
            n = t // (t & -t)
            if n < m:
                res = ("D", s); break
            if n == m:
                res = ("F", s); break
        out[m] = res
    return out


def read_dump(path):
    d = {}
    with open(path) as f:
        for line in f:
            m, o, s = line.split()
            d[int(m)] = (o, int(s))
    return d


def main():
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8")

    def log(s):
        print(s, flush=True); logf.write(s + "\n")

    started = datetime.now(timezone.utc).isoformat()
    log(f"=== B3 run {started}  python {sys.version.split()[0]}")
    summary = {"started": started, "tests": []}
    ok_all = True

    cases = [(5, 100, None), (3, 60, None), (181, 10, None), (7, 2963, None),
             (15, 2437, None), (31, 3645, None), (7, 32603, 8_000)]
    expected_found = {5: {13, 17}, 181: {27, 35}, 3: set()}
    for q, C, mhi in cases:
        t0 = time.time()
        dump = os.path.join(OUT, f"dump_q{q}_C{C}.txt")
        line = run_c(q, C, mhi, dump)
        cd = read_dump(dump)
        t1 = time.time()
        dump_ref = os.path.join(OUT, f"dump_q{q}_C{C}_noblock.txt")
        run_c(q, C, mhi, dump_ref, noblock=True)
        block_equal = read_dump(dump_ref) == cd
        t1b = time.time()
        pr = python_reference(q, C, mhi)
        t2 = time.time()
        same_keys = set(cd) == set(pr)
        mism = [m for m in pr if cd.get(m) != pr[m]]
        found = {m for m, (o, _) in cd.items() if o == "F"}
        ok = same_keys and not mism and block_equal
        if q in expected_found:
            ok = ok and found == expected_found[q]
        if (q, C) == (7, 2963):
            ok = ok and len(cd) == 74277
        ok_all &= ok
        rec = {"q": q, "C": C, "mhi": mhi, "candidates": len(cd), "found": sorted(found),
               "candidate_sets_equal": same_keys, "mismatches": len(mism),
               "block_equals_noblock": block_equal,
               "c_seconds": round(t1 - t0, 2), "c_noblock_seconds": round(t1b - t1, 2),
               "python_seconds": round(t2 - t1b, 2), "pass": ok}
        summary["tests"].append(rec)
        log(("PASS " if ok else "FAIL ") + json.dumps(rec))
        if mism:
            log(f"   first mismatches: {[(m, cd.get(m), pr[m]) for m in mism[:5]]}")
    summary["all_pass"] = ok_all
    summary["finished"] = datetime.now(timezone.utc).isoformat()
    json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    log(f"=== ALL PASS: {ok_all}")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
