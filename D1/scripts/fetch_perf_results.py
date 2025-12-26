import os
import re
import csv
import numpy as np

RESULTS_DIR = "results"
OUTPUT_CSV = os.path.join(RESULTS_DIR, "perf_summary.csv")

MATRIX_DIRS = {"fd12", "sinc12", "c-46", "epb1", "human_gene2"}
SCHEDULINGS = ["static", "dynamic", "guided"]
MEM_SIZES = ["32GB", "64GB", "128GB"]
THREAD_COUNTS = [1, 2, 4, 8, 16, 32, 64]

def parse_perf_file(filepath):
    executions = []
    with open(filepath, "r") as f:
        lines = [line.rstrip() for line in f]
    chunks = []
    current = []
    for line in lines:
        if line.startswith("# Execution:"):
            if current:
                chunks.append(current)
                current = []
        elif line.strip() == "# END EXECUTION":
            if current:
                chunks.append(current)
                current = []
        else:
            current.append(line)
    if current:
        chunks.append(current)

    def extract_value(pattern, chunk):
        for line in chunk:
            if pattern in line:
                parts = re.split(r"\s+", line.strip())
                try:
                    value = int(parts[0].replace(",", ""))
                except ValueError:
                    value = 0
                return value
        return 0

    for chunk in chunks:
        cache_refs = extract_value("cache-references", chunk)
        cache_misses = extract_value("cache-misses", chunk)
        l1_loads = extract_value("L1-dcache-loads", chunk)
        l1_misses = extract_value("L1-dcache-load-misses", chunk)
        llc_loads = extract_value("LLC-loads", chunk)
        llc_misses = extract_value("LLC-load-misses", chunk)

        if cache_refs and l1_loads and llc_loads:
            results = {
                "cache_misses_percent": 100 * cache_misses / cache_refs if cache_refs else np.nan,
                "l1_misses_percent": 100 * l1_misses / l1_loads if l1_loads else np.nan,
                "llc_misses_percent": 100 * llc_misses / llc_loads if llc_loads else np.nan,
            }
            executions.append(results)

    return executions


def main():
    rows = []
    for matrix in sorted(os.listdir(RESULTS_DIR)):
        matrix_path = os.path.join(RESULTS_DIR, matrix)
        if not os.path.isdir(matrix_path) or matrix not in MATRIX_DIRS:
            continue

        for sched in SCHEDULINGS:
            sched_path = os.path.join(matrix_path, sched)
            if not os.path.isdir(sched_path):
                continue

            for mem in MEM_SIZES:
                perf_path = os.path.join(sched_path, mem, "perf")
                if not os.path.isdir(perf_path):
                    continue

                for t in THREAD_COUNTS:
                    pattern = re.compile(rf"{matrix}_perf_NT{t}\.txt$")
                    matches = [f for f in os.listdir(perf_path) if pattern.match(f)]
                    if not matches:
                        continue
                    filepath = os.path.join(perf_path, matches[0])

                    executions = parse_perf_file(filepath)
                    if not executions:
                        continue

                    # Average percentages over runs
                    avg_cache = np.mean([e["cache_misses_percent"] for e in executions])
                    avg_l1 = np.mean([e["l1_misses_percent"] for e in executions])
                    avg_llc = np.mean([e["llc_misses_percent"] for e in executions])

                    rows.append({
                        "matrix": matrix,
                        "scheduling": sched,
                        "memory": mem,
                        "threads": t,
                        "cache_misses_percent": round(avg_cache, 3),
                        "l1_misses_percent": round(avg_l1, 3),
                        "llc_misses_percent": round(avg_llc, 3),
                    })

    # Write CSV
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as csvfile:
        fieldnames = [
            "matrix",
            "scheduling",
            "memory",
            "threads",
            "cache_misses_percent",
            "l1_misses_percent",
            "llc_misses_percent",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (r["matrix"], r["scheduling"], r["memory"], r["threads"])))

    print(f"Perf summary CSV written to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
