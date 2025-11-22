import os
import re
import csv
import numpy as np

RESULTS_DIR = "results/cachegrind_logs"
OUTPUT_CSV = os.path.join("results", "cachegrind_summary.csv")

MATRIX_DIRS = {"fd12", "sinc12", "c-46", "epb1", "human_gene2"}
SCHEDULINGS = ["static", "dynamic", "guided"]

PATTERNS = {
    "I_refs": re.compile(r"I\s+refs:\s+([\d,]+)"),
    "I1_misses": re.compile(r"I1\s+misses:\s+([\d,]+)"),
    "LLi_misses": re.compile(r"LLi\s+misses:\s+([\d,]+)"),
    "D_refs": re.compile(r"D\s+refs:\s+([\d,]+)"),
    "D1_misses": re.compile(r"D1\s+misses:\s+([\d,]+)"),
    "LLd_misses": re.compile(r"LLd\s+misses:\s+([\d,]+)"),
    "LL_refs": re.compile(r"LL\s+refs:\s+([\d,]+)"),
    "LL_misses": re.compile(r"LL\s+misses:\s+([\d,]+)"),
}

def parse_err_file(filepath):
    counts = {}
    with open(filepath, "r") as f:
        content = f.read()

    for key, pattern in PATTERNS.items():
        match = pattern.search(content)
        if match:
            counts[key] = int(match.group(1).replace(",", ""))
        else:
            counts[key] = 0

    metrics = {}
    metrics["I1_miss_percent"] = (counts["I1_misses"] / counts["I_refs"] * 100) if counts["I_refs"] else 0.0
    metrics["LLi_miss_percent"] = (counts["LLi_misses"] / counts["I_refs"] * 100) if counts["I_refs"] else 0.0
    metrics["D1_miss_percent"] = (counts["D1_misses"] / counts["D_refs"] * 100) if counts["D_refs"] else 0.0
    metrics["LLd_miss_percent"] = (counts["LLd_misses"] / counts["D_refs"] * 100) if counts["D_refs"] else 0.0
    metrics["LL_miss_percent"] = (counts["LL_misses"] / counts["LL_refs"] * 100) if counts["LL_refs"] else 0.0

    return metrics

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

            err_files = [f for f in os.listdir(sched_path) if f.endswith(".err")]
            if not err_files:
                continue

            err_path = os.path.join(sched_path, err_files[0])
            metrics = parse_err_file(err_path)

            row = {"matrix": matrix, "scheduling": sched}
            row.update({k: round(v, 3) for k, v in metrics.items()})
            rows.append(row)

    # Write CSV
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as csvfile:
        fieldnames = ["matrix", "scheduling",
                      "I1_miss_percent", "LLi_miss_percent",
                      "D1_miss_percent", "LLd_miss_percent", "LL_miss_percent"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (r["matrix"], r["scheduling"])))

    print(f"Cachegrind summary CSV written to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
