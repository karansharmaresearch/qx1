import bisect, json, os, subprocess, sys, tempfile, time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "code"))
import bounds as b

EXE = os.path.join(ROOT, "code", "cycles.exe")
PROD = [(7, 91312), (15, 44001), (31, 31765), (63, 12631), (127, 67960),
        (255, 94924), (511, 83317), (1023, 94354)]


def admissible(q, lo, hi):
    adm, x = set(), 1
    for _ in range(q):
        adm.add(x); x = 2 * x % q
    return [m for m in range(lo | 1, hi, 2) if m % q in adm]


def windows(q, B, nc):
    span = 2 * q * nc + 2 * q
    out = []
    for lo in (3, max(3, B // 2), max(3, B - span)):
        c = admissible(q, lo, min(B, lo + span))[:nc]
        out.append((c[0], c[-1] + 1))
    return out


def python_ref(q, C, ms, BA):
    pairs = sorted((BA[c], c) for c in range(1, C + 1))
    keys = [p[0] for p in pairs]
    suf = [0] * (len(pairs) + 1)
    for i in range(len(pairs) - 1, -1, -1):
        suf[i] = max(suf[i + 1], pairs[i][1])
    out = {}
    for m in ms:
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


def c_dump(q, C, lo, hi, noblock):
    st = os.path.join(ROOT, "data", "bounds", f"staircase_q{q}_C{C}.txt")
    with tempfile.TemporaryDirectory() as td:
        dump = os.path.join(td, "d.txt")
        args = [EXE, st, td, "--mlo", str(lo), "--mhi", str(hi), "--dump", dump] + (["--noblock"] if noblock else [])
        r = subprocess.run(args, capture_output=True, text=True)
        if r.returncode not in (0, 3):
            raise RuntimeError(r.stderr)
        d = {}
        for line in open(dump):
            m, o, s = line.split()
            d[int(m)] = (o, int(s))
    return d


def main():
    nc = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    report, ok_all = [], True
    for q, C in PROD:
        t0 = time.time()
        _, BA = b.bounds_exact(q, C)
        B = max(BA[1:])
        for lo, hi in windows(q, B, nc):
            blk, nob = c_dump(q, C, lo, hi, False), c_dump(q, C, lo, hi, True)
            py = python_ref(q, C, admissible(q, lo, hi), BA)
            ok = blk == nob == py and len(py) > 0 and all(o != "F" for o, _ in py.values())
            ok_all &= ok
            rec = {"q": q, "C": C, "mlo": lo, "mhi": hi, "candidates": len(py),
                   "outcomes": {k: sum(1 for o, _ in py.values() if o == k) for k in "DSF"},
                   "max_steps": max(s for _, s in py.values()),
                   "block_eq_noblock": blk == nob, "c_eq_python": blk == py, "pass": ok}
            report.append(rec)
            print(("PASS " if ok else "FAIL ") + json.dumps(rec), flush=True)
        print(f"   q={q} done in {time.time() - t0:.0f}s", flush=True)
    out = os.path.join(ROOT, "data", "checks", "production_crosscheck.json")
    json.dump({"nc": nc, "cases": report, "all_pass": ok_all,
               "exe_note": "code/cycles.exe (the production binary)"}, open(out, "w"), indent=1)
    with open(os.path.join(ROOT, "data", "checks", "log.txt"), "a") as f:
        f.write(f"=== production cross-check nc={nc}: all_pass={ok_all}\n")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
