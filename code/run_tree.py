import hashlib, json, os, platform, subprocess, sys
from datetime import datetime, timezone

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXE = os.path.join(ROOT, "code", "tree.exe")
OUT = os.path.join(ROOT, "data", "tree")


def main():
    q, D = int(sys.argv[1]), int(sys.argv[2])
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"q{q}_D{D}.json")
    started = datetime.now(timezone.utc).isoformat()
    p = subprocess.run([EXE, str(q), str(D), path], capture_output=True, text=True)
    if p.returncode:
        sys.exit(p.stderr)
    d = json.load(open(path))
    d["provenance"] = {
        "started": started, "finished": datetime.now(timezone.utc).isoformat(),
        "exe_sha256": hashlib.sha256(open(EXE, "rb").read()).hexdigest(),
        "machine": {"platform": platform.platform(), "processor": platform.processor(),
                    "cpu_count": os.cpu_count()},
    }
    json.dump(d, open(path, "w"), indent=1)
    with open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8") as f:
        f.write(f"{started} {p.stdout.strip()} exe={d['provenance']['exe_sha256'][:12]}\n")
    print(p.stdout.strip())


if __name__ == "__main__":
    main()
