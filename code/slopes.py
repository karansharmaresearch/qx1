import glob, json, math, os, re

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "data", "slopes")
WIN = {7: 4, 15: 6, 31: 10, 63: 15}
MARGIN = {7: 10, 15: 10, 31: 20, 63: 20}


def gamma(q):
    lo, hi = 1e-12, 0.999999
    for _ in range(200):
        g = (lo + hi) / 2
        if 1 + q ** (g - 1) - 2 ** g > 0:
            lo = g
        else:
            hi = g
    return (lo + hi) / 2


def N(H, a, b):
    return sum(H[i][j] for i in range(a + 1) for j in range(b + 1))


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {}
    for q in sorted(WIN):
        files = glob.glob(os.path.join(ROOT, "data", "tree", f"q{q}_D*.json"))
        D, path = max((int(re.search(r"_D(\d+)", f).group(1)), f) for f in files)
        d = json.load(open(path))
        H, w = d["H"], WIN[q]
        amax = D - MARGIN[q]
        wins = []
        for a in range(w, amax - w + 1, w):
            n1, n2 = N(H, a, D), N(H, a + w, D)
            if n1 > 0:
                wins.append((a, a + w, math.log(n2 / n1) / (w * math.log(10))))
        last = [s for _, _, s in wins[-3:]]
        up = [x for x in wins if x[0] >= amax / 2] or wins[-1:]
        ups = [s for _, _, s in up]
        ladder = {a: [N(H, a, b) if b <= D else None for b in (a, a + 10, a + 20)] + [N(H, a, D)]
                  for a in range(w, amax + 1, w)}
        stab = next(m for m in range(0, MARGIN[q] + 1)
                    if all(N(H, a, b) == N(H, a, D) for a in range(1, amax + 1) for b in range(a + m, D + 1)))
        res[q] = {"D": D, "margin": MARGIN[q], "observed_stability_margin": stab, "window": w,
                  "nodes": d["nodes"], "gamma": gamma(q), "windows": wins,
                  "last3": {"mean": sum(last) / len(last), "halfrange": (max(last) - min(last)) / 2},
                  "upper": {"windows": up, "min": min(ups), "max": max(ups),
                            "mean": sum(ups) / len(ups), "halfrange": (max(ups) - min(ups)) / 2},
                  "ladder": ladder, "source": os.path.relpath(path, ROOT).replace("\\", "/")}
    json.dump(res, open(os.path.join(OUT, "fits.json"), "w"), indent=1)

    for q, r in res.items():
        print(f"q={q:3d} D={r['D']} upper-half slopes {r['upper']['min']:.4f}-{r['upper']['max']:.4f} "
              f"gamma={r['gamma']:.4f} observed stability margin={r['observed_stability_margin']}")


if __name__ == "__main__":
    main()
