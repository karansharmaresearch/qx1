import hashlib, json, os, platform, subprocess, sys
from datetime import datetime, timezone

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXE = os.path.join(ROOT, "code", "cycles.exe")


def main():
    q, C = int(sys.argv[1]), int(sys.argv[2])
    threads = sys.argv[3] if len(sys.argv) > 3 else str(os.cpu_count())
    st = os.path.join(ROOT, "data", "bounds", f"staircase_q{q}_C{C}.txt")
    if not os.path.exists(st):
        subprocess.run([sys.executable, os.path.join(ROOT, "code", "make_staircase.py"),
                        str(q), str(C)], check=True)
    out = os.path.join(ROOT, "data", "cycles", f"q{q}_C{C}")
    os.makedirs(out, exist_ok=True)
    body = "".join(l for l in open(st) if not l.startswith("#"))
    meta = {
        "q": q, "C": C, "threads": int(threads),
        "staircase": os.path.relpath(st, ROOT).replace("\\", "/"),
        "staircase_body_sha256": hashlib.sha256(body.encode()).hexdigest(),
        "exe_sha256": hashlib.sha256(open(EXE, "rb").read()).hexdigest(),
        "machine": {"platform": platform.platform(), "processor": platform.processor(),
                    "cpu_count": os.cpu_count()},
        "started": datetime.now(timezone.utc).isoformat(),
    }
    logpath = os.path.join(out, "run.log")
    previous = []
    if os.path.exists(logpath):
        for line in open(logpath, encoding="utf-8"):
            if line.startswith("=== start "):
                previous.append(json.loads(line.split(" ", 3)[3]))
    for prev in previous:
        for key in ("exe_sha256", "staircase_body_sha256"):
            if prev.get(key) != meta[key]:
                sys.exit(f"refusing to resume: {key} differs from an earlier invocation "
                         f"({prev.get(key)} vs {meta[key]})")
    meta["previous_invocations"] = [p["started"] for p in previous]
    log = open(logpath, "a", encoding="utf-8")
    log.write(f"=== start {meta['started']} {json.dumps(meta)}\n"); log.flush()
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    p = subprocess.run([EXE, st.replace("\\", "/"), out.replace("\\", "/")],
                       env=env, capture_output=True, text=True)
    log.write(p.stderr); log.write(p.stdout)
    meta["finished"] = datetime.now(timezone.utc).isoformat()
    meta["exit_code"] = p.returncode
    meta["summary"] = p.stdout.strip()
    res = [f for f in os.listdir(out) if f.startswith("result_")]
    if res:
        txt = open(os.path.join(out, res[0])).read()
        try:
            meta["result"] = json.loads(txt)
        except json.JSONDecodeError:
            meta["result"] = json.loads(txt.replace("\\", "/"))
    json.dump(meta, open(os.path.join(out, "run.json"), "w"), indent=1)
    log.write(f"=== end {meta['finished']} exit={p.returncode}\n")
    print(meta["summary"]); print("exit", p.returncode)
    sys.exit(0 if p.returncode in (0, 3) else p.returncode)


if __name__ == "__main__":
    main()
