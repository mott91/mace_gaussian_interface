"""Summarise the link-717 core-count test (thesis docs/handoff/l717-cores-prompt.md).

usage: python analyse.py <test_dir> <campaign_log>

For each <test_dir>/npN run: link-717 cpu and elap from the #P "Leave Link" line, whole-job
cpu and elapsed, the thread samples of l717.exe (threads, threads above 50 % CPU in top's
instantaneous frame), and whether every "Fundamental Bands" table matches the campaign log
character for character.
"""

import re
import sys
from pathlib import Path

LEAVE = re.compile(r"Leave Link\s+(\d+) at .*cpu:\s+([\d.]+) elap:\s+([\d.]+)")
JOB = re.compile(r"(Job cpu time|Elapsed time):\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+([\d.]+)")


def seconds(d, h, m, s):
    return ((int(d) * 24 + int(h)) * 60 + int(m)) * 60 + float(s)


def fundamentals(text):
    """All 'Fundamental Bands' tables, each as a tuple of stripped rows."""
    lines = text.splitlines()
    tables = []
    for i, line in enumerate(lines):
        if line.strip().endswith("Fundamental Bands") and lines[i + 1].strip().startswith("---"):
            rows = []
            for row in lines[i + 3 :]:
                if not re.match(r"^\s+\d+\(\d+\)", row):
                    break
                rows.append(row.strip())
            tables.append(tuple(rows))
    return tables


def samples(path):
    if not path.exists():
        return None
    blocks = path.read_text().split("### ")[1:]
    threads, busy, pcpu = [], [], []
    for b in blocks:
        m = re.search(r"Threads:\s+(\d+)", b)
        if m:
            threads.append(int(m.group(1)))
        top = b.split("instantaneous)")[-1]
        cpus = [float(x) for x in re.findall(r"^\s*\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S\s+([\d.]+)", top, re.M)]
        busy.append(sum(c > 50 for c in cpus))
        pcpu.append(sum(cpus))
    aff = sorted(set(re.findall(r"Cpus_allowed_list:\s+(\S+)", path.read_text())))
    omp = sorted(set(re.findall(r"OMP_NUM_THREADS=(\d+)", path.read_text())))
    return dict(n=len(blocks), threads=sorted(set(threads)), busy=sorted(set(busy)),
                pcpu=(min(pcpu), max(pcpu)) if pcpu else None, aff=aff, omp=omp)


def main():
    test_dir, campaign_log = Path(sys.argv[1]), Path(sys.argv[2])
    ref = fundamentals(campaign_log.read_text(errors="replace"))
    print(f"campaign: {len(ref)} Fundamental Bands tables, {len(ref[0])} rows each\n")
    print("| N | l717 elap / s | l717 cpu / s | ratio | job elap / s | job cpu / s | ratio | "
          "l717 threads | busy (>50 %) | sum %CPU | OMP_NUM_THREADS | affinity | samples | "
          "fundamentals = campaign |")
    print("|" + "---|" * 14)
    runs = sorted(test_dir.glob("np*/p.log"), key=lambda p: int(p.parent.name[2:]))
    for log in runs:
        text = log.read_text(errors="replace")
        l717 = [(float(c), float(e)) for k, c, e in LEAVE.findall(text) if k == "717"]
        job = {k: seconds(*v) for k, *v in JOB.findall(text)}
        s = samples(log.parent / "samples.txt")
        same = fundamentals(text) == ref
        c, e = l717[-1] if l717 else (float("nan"), float("nan"))
        jc, je = job.get("Job cpu time", float("nan")), job.get("Elapsed time", float("nan"))
        pc = f"{s['pcpu'][0]:.0f}–{s['pcpu'][1]:.0f}" if s and s["pcpu"] else "–"
        print(f"| {log.parent.name[2:]} | {e:.1f} | {c:.1f} | {c / e:.2f} | {je:.1f} | {jc:.1f} | "
              f"{jc / je:.2f} | {s['threads'] if s else '–'} | {s['busy'] if s else '–'} | {pc} | "
              f"{','.join(s['omp']) if s else '–'} | {','.join(s['aff']) if s else '–'} | "
              f"{s['n'] if s else 0} | {'yes' if same else 'NO'} |")


if __name__ == "__main__":
    main()
