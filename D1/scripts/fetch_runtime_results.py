import os
import re
import csv
import numpy as np

# Results directory
RESULTS_DIR = "results"
OUTPUT_CSV = os.path.join(RESULTS_DIR, "runtime_summary.csv")

# Matrix folders
MATRIX_DIRS = {"fd12", "sinc12", "c-46", "epb1", "human_gene2"}

SCHEDULINGS = ["static", "dynamic", "guided"]
MEM_SIZES = ["32GB", "64GB", "128GB"]
THREAD_COUNTS = [1, 2, 4, 8, 16, 32, 64]

def extract_times(filepath):
    times = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if re.fullmatch(r"\d+", line):
                continue
            try:
                times.append(float(line))
            except ValueError:
                pass
    return times

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
                mem_path = os.path.join(sched_path, mem, "runtime")
                if not os.path.isdir(mem_path):
                    continue

                for t in THREAD_COUNTS:
                    filename_pattern = re.compile(rf"{matrix}_NT{t}\.txt$")
                    # Find the correct runtime file
                    matches = [f for f in os.listdir(mem_path) if filename_pattern.match(f)]
                    if not matches:
                        continue
                    filepath = os.path.join(mem_path, matches[0])

                    times = extract_times(filepath)
                    if not times:
                        continue
                    p90 = np.percentile(times, 90)

                    rows.append({
                        "matrix": matrix,
                        "scheduling": sched,
                        "memory": mem,
                        "threads": t,
                        "p90_runtime_ms": round(p90, 6)
                    })

    # Write results
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as csvfile:
        fieldnames = ["matrix", "scheduling", "memory", "threads", "p90_runtime_ms"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (r["matrix"], r["scheduling"], r["memory"], r["threads"])))

    print(f"CSV file written to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
