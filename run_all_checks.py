"""run_all_checks.py -- run every check script in scripts/ and report.

Usage:
    python3 run_all_checks.py            # one script at a time
    python3 run_all_checks.py -j 4       # four at a time (needs about 2 GB of memory per job for the largest scripts)
    python3 run_all_checks.py NAME ...   # only the scripts whose file name contains one of the given words

Each script runs in its own process, from inside scripts/. A script is CLEAN when it exits with code 0,
prints no line starting with FAIL and prints no Traceback. The lines it prints that start with PASS or FAIL
are counted. A report is written to CHECK_REPORT.md. The exit code is 0 only if every script is clean.

Needs Python 3 with numpy, scipy, sympy, mpmath and matplotlib (tested with Python 3.13).
"""
import os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "scripts")
TIMEOUT = 1800  # seconds for one script


def run(name):
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, name], cwd=SCRIPTS, capture_output=True, text=True, timeout=TIMEOUT)
        out, code = p.stdout + p.stderr, p.returncode
    except subprocess.TimeoutExpired:
        out, code = "TIMEOUT", -1
    lines = [l.lstrip() for l in out.splitlines()]
    npass = sum(l.startswith("PASS") for l in lines)
    nfail = sum(l.startswith("FAIL") for l in lines)
    clean = code == 0 and nfail == 0 and "Traceback" not in out
    return dict(name=name, code=code, sec=time.time() - t0, npass=npass, nfail=nfail, clean=clean, tail=out[-400:])


def main():
    args = sys.argv[1:]
    jobs = 1
    if "-j" in args:
        i = args.index("-j")
        jobs = int(args[i + 1])
        del args[i:i + 2]
    names = sorted(f for f in os.listdir(SCRIPTS) if f.endswith(".py"))
    if args:
        names = [n for n in names if any(a in n for a in args)]
    print(f"running {len(names)} scripts from scripts/, {jobs} at a time")
    results = []
    with ThreadPoolExecutor(jobs) as ex:
        for r in ex.map(run, names):
            results.append(r)
            mark = "clean " if r["clean"] else "NOT CLEAN"
            print(f"{mark:9s} {r['npass']:3d} pass {r['nfail']:2d} fail {r['sec']:7.1f} s  {r['name']}", flush=True)
    bad = [r for r in results if not r["clean"]]
    quiet = [r for r in results if r["clean"] and r["npass"] == 0]
    total_pass = sum(r["npass"] for r in results)
    total_fail = sum(r["nfail"] for r in results)
    summary = (f"{len(results)} scripts: {len(results) - len(bad)} clean, {len(bad)} not clean; "
               f"{total_pass} PASS lines, {total_fail} FAIL lines; "
               f"{len(quiet)} clean scripts print numbers or figures without PASS lines")
    print("\n" + summary)
    with open(os.path.join(HERE, "CHECK_REPORT.md"), "w", encoding="utf-8") as f:
        f.write("# Check report\n\n" + summary + "\n\n")
        f.write(f"Python {sys.version.split()[0]}, run on {time.strftime('%Y-%m-%d')}.\n\n")
        f.write("| script | result | PASS | FAIL | seconds |\n|---|---|---|---|---|\n")
        for r in results:
            f.write(f"| {r['name']} | {'clean' if r['clean'] else 'NOT CLEAN'} | {r['npass']} | {r['nfail']} | {r['sec']:.1f} |\n")
        if bad:
            f.write("\n## Not clean: the end of each output\n\n")
            for r in bad:
                f.write(f"### {r['name']} (exit code {r['code']})\n\n```\n{r['tail']}\n```\n\n")
    for r in bad:
        print("NOT CLEAN:", r["name"], "exit code", r["code"])
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
