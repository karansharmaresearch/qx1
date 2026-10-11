import hashlib, json, os, sys

sys.path.insert(0, os.path.dirname(__file__))
import bounds as b

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def main():
    q, C = int(sys.argv[1]), int(sys.argv[2])
    outdir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "data", "bounds")
    J, Bup, gap = b.certified_bounds(q, C)
    st = b.staircase(Bup, C)
    sub, x = {1}, 2 % q
    while x != 1:
        sub.add(x)
        x = 2 * x % q
    sub = sorted(sub)
    body = "".join(f"{c} {B}\n" for c, B in st)
    path = os.path.join(outdir, f"staircase_q{q}_C{C}.txt")
    with open(path, "w", newline="\n") as f:
        f.write(f"# q {q}\n# C {C}\n# B {st[0][1]}\n")
        f.write("# subgroup " + " ".join(map(str, sub)) + "\n")
        f.write(f"# entries {len(st)}\n")
        f.write(body)
    meta = {"q": q, "C": C, "B": st[0][1], "entries": len(st), "subgroup": sub,
            "sha256": hashlib.sha256(body.encode()).hexdigest()}
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
