import json, os, platform, sys, time
from datetime import datetime, timezone

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def dec(x):
    return 0 if x <= 1 else len(str(x - 1))


def T(q, n):
    t = q * n + 1
    return t // (t & -t)


def run(q, D):
    k = q.bit_length()
    assert q == (1 << k) - 1, "q must be 2^k - 1"
    L = 10 ** D
    lim = q * L + 1
    dlog = {}
    x = 1
    for a in range(k):
        dlog[x] = a
        x = 2 * x % q
    H = [[0] * (D + 1) for _ in range(D + 1)]
    H[0][0] = 1
    nodes, checked = 1, 0
    stack = [(1, 1)]
    while stack:
        m, Mm = stack.pop()
        a = dlog.get(m % q)
        if a is None:
            continue
        v = (k - a % k) % k or k
        cur = m << v
        while cur <= lim:
            n = (cur - 1) // q
            if n != 1:
                Mn = n if n > Mm else Mm
                y, My = n, n
                while y != 1:
                    y = T(q, y)
                    if y > My:
                        My = y
                assert My == Mn and Mn <= L, (q, n, My, Mn)
                checked += 1
                H[dec(n)][dec(Mn)] += 1
                nodes += 1
                stack.append((n, Mn))
            cur <<= k
    return H, nodes, checked


def main():
    q, D = int(sys.argv[1]), int(sys.argv[2])
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "data", "tree", f"q{q}_D{D}.json")
    t0 = time.time()
    H, nodes, checked = run(q, D)
    sec = time.time() - t0
    d = {"q": q, "k": q.bit_length(), "D": D, "Lmax": str(10 ** D), "nodes": nodes,
         "forward_certified": checked + 1, "seconds": round(sec, 3), "implementation": "bigtree.py",
         "hist_note": "H[a][b] = #{n: dec(n)=a, dec(orbit max)=b}, dec(x)=min{a: x<=10^a}",
         "H": H,
         "provenance": {"finished": datetime.now(timezone.utc).isoformat(),
                        "python": sys.version.split()[0], "platform": platform.platform()}}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(d, open(out, "w"), indent=1)
    print(f"q={q} D={D} nodes={nodes} certified={checked + 1} {sec:.1f}s")


if __name__ == "__main__":
    main()
