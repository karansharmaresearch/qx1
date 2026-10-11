import os
from mpmath import mp, iv, mpf, floor

N = 100
mp.dps = N + 40
iv.dps = N + 40


def f(q, x):
    return 1 + iv.mpf(q) ** (x - 1) - iv.mpf(2) ** x


def gamma(q):
    lo, hi = mpf("0.001"), mpf("0.999")
    for _ in range(4 * (N + 30)):
        mid = (lo + hi) / 2
        if 1 + mpf(q) ** (mid - 1) - mpf(2) ** mid > 0:
            lo = mid
        else:
            hi = mid
    d = int(floor(lo * 10 ** N))
    a = iv.mpf(d) / 10 ** N
    b = iv.mpf(d + 1) / 10 ** N
    assert f(q, a).a > 0 and f(q, b).b < 0
    return "0." + str(d).zfill(N)


out = os.path.join(os.path.dirname(__file__), "..", "data", "gamma")
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "gamma.txt"), "w") as g:
    for q in [7, 15, 31, 63, 127, 255, 511, 1023]:
        s = gamma(q)
        print(q, s)
        g.write(f"{q} {s}\n")
