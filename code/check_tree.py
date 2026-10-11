import json, os, subprocess, sys, tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXE = os.path.join(ROOT, "code", "tree.exe")
OUT = os.path.join(ROOT, "data", "tree")


def T(q, n):
    t = q * n + 1
    return t // (t & -t)


def dec(x):
    a, p = 0, 1
    while p < x:
        a += 1; p *= 10
    return a


def python_hist(q, D):
    L = 10 ** D
    seen, frontier = {1}, [1]
    while frontier:
        nxt = []
        for m in frontier:
            v = 1
            while (2 ** v * m - 1) // q <= L:
                t = 2 ** v * m - 1
                if t % q == 0:
                    n = t // q
                    if n not in seen:
                        seen.add(n); nxt.append(n)
                v += 1
        frontier = nxt
    H = [[0] * (D + 1) for _ in range(D + 1)]
    for n in seen:
        M, x = n, n
        while x != 1:
            x = T(q, x)
            if x > M:
                M = x
        if M > L:
            raise AssertionError(f"q={q}: tree member {n} has orbit max {M} > L")
        H[dec(n)][dec(M)] += 1
    return H, len(seen)


def main():
    os.makedirs(OUT, exist_ok=True)
    cases = [(7, 12), (7, 16), (15, 30), (31, 36), (63, 36)]
    report, ok_all = [], True
    for q, D in cases:
        with tempfile.TemporaryDirectory() as td:
            out = os.path.join(td, "t.json")
            subprocess.run([EXE, str(q), str(D), out], check=True, capture_output=True)
            c = json.load(open(out))
        H, n = python_hist(q, D)
        ok = (H == c["H"]) and (n == c["nodes"])
        ok_all &= ok
        rec = {"q": q, "D": D, "nodes_c": c["nodes"], "nodes_py": n, "hist_equal": H == c["H"], "pass": ok}
        report.append(rec)
        print(("PASS " if ok else "FAIL ") + json.dumps(rec), flush=True)
    json.dump({"cases": report, "all_pass": ok_all}, open(os.path.join(OUT, "verify.json"), "w"), indent=1)
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
