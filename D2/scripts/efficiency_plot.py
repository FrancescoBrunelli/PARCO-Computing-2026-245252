import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RESULTS_DIR = "results"
PLOTS_DIR = "plots"
os.makedirs(PLOTS_DIR, exist_ok=True)

# Partitioning types and colors
PARTITIONS = {
    "1D_Block": "blue",
    "1D_Cyclic": "red",
    "2D_Block": "green"
}

# MPI process counts
PROCESSES = [1, 2, 4, 8, 16, 32, 64, 128]

def extract_csv_p90(filepath):
    df = pd.read_csv(filepath, comment="#", header=None)
    times = df[0].values   # Execution_time[ms]
    return np.percentile(times, 90)

def collect_matrix(matrix_name):
    matrix_path = os.path.join(RESULTS_DIR, matrix_name)

    # Sequential baseline
    seq_file = os.path.join(matrix_path, f"{matrix_name}_NP1.csv")
    if not os.path.isfile(seq_file):
        print(f"ERROR: missing baseline {seq_file}")
        return None

    T1 = extract_csv_p90(seq_file)

    data = {part: {} for part in PARTITIONS}

    for part in PARTITIONS:
        part_path = os.path.join(matrix_path, part)
        if not os.path.isdir(part_path):
            continue

        for p in PROCESSES[1:]:
            csv_file = os.path.join(part_path, f"{matrix_name}_NP{p}.csv")
            if not os.path.isfile(csv_file):
                print(f"Warning: missing {csv_file}")
                continue

            Tp = extract_csv_p90(csv_file)
            efficiency = T1 / (p * Tp)
            data[part][p] = efficiency

    return data

def plot_efficiency(matrix_name, data):
    plt.figure(figsize=(8, 5))

    # Equally spaced x-axis
    x_indices = np.arange(len(PROCESSES) - 1)   # skip NP=1

    for part, color in PARTITIONS.items():
        y = [data[part].get(p, np.nan) for p in PROCESSES[1:]]
        plt.plot(x_indices, y, marker="o", color=color, label=part)

    plt.xticks(x_indices, PROCESSES[1:])
    plt.xlabel("Number of MPI Processes")
    plt.ylabel("Parallel Efficiency")
    plt.title(f"{matrix_name} - Efficiency")
    plt.ylim(bottom=0)

    # Horizontal reference lines
    ax = plt.gca()
    for y in ax.get_yticks():
        plt.axhline(y=y, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

    plt.legend()
    plt.tight_layout()

    out = os.path.join(PLOTS_DIR, f"{matrix_name}_efficiency.png")
    plt.savefig(out)
    plt.close()

    print(f"Saved {out}")

def main():
    for matrix in os.listdir(RESULTS_DIR):
        matrix_path = os.path.join(RESULTS_DIR, matrix)
        if not os.path.isdir(matrix_path):
            continue

        # Only matrices with NP1 baseline
        baseline = os.path.join(matrix_path, f"{matrix}_NP1.csv")
        if not os.path.isfile(baseline):
            continue

        print(f"Processing {matrix}")
        data = collect_matrix(matrix)
        if data is None:
            continue

        plot_efficiency(matrix, data)

if __name__ == "__main__":
    main()
