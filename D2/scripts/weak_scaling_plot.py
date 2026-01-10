import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Only plot this matrix
MATRIX_NAME = "weak_scaling"

RESULTS_DIR = "results"
PLOTS_DIR = "plots"
os.makedirs(PLOTS_DIR, exist_ok=True)

# Partitioning types and colors
PARTITIONS = {
    "1D_Block": "blue",
    "1D_Cyclic": "red",
    "2D_Block": "green"
}

# Metrics:
# (column_index, title_label (no units), y_label (with units), filename_suffix)
METRICS = [
    (0, "Execution Time", "Execution Time [ms]", "execution_time"),
    (4, "GFLOPS", "GFLOPS", "gflops"),
    (5, "Communication Time", "Communication Time [s]", "communication_time"),
    (6, "Memory Usage", "Memory Usage [MB]", "memory_usage")
]

# MPI processes
PROCESSES = [2, 4, 8, 16, 32, 64, 128]

def extract_csv_p90(filepath, col_idx, comm_time=False):
    df = pd.read_csv(filepath, comment="#", header=None)
    values = df[col_idx].values
    p90 = np.percentile(values, 90)
    if comm_time:
        p90 /= 1000.0   # ms → s
    return p90

def collect_data(matrix_path):
    data = {col_idx: {} for col_idx, _, _, _ in METRICS}

    for col_idx, _, _, _ in METRICS:
        for part in PARTITIONS:
            data[col_idx][part] = {}

    for partition in PARTITIONS:
        part_path = os.path.join(matrix_path, partition)
        if not os.path.isdir(part_path):
            continue

        for p in PROCESSES:
            csv_file = os.path.join(part_path, f"{MATRIX_NAME}_NP{p}.csv")
            if not os.path.isfile(csv_file):
                print(f"Warning: missing {csv_file}")
                continue

            for col_idx, _, _, _ in METRICS:
                try:
                    val = extract_csv_p90(csv_file, col_idx, col_idx == 5)
                    data[col_idx][partition][p] = val
                except Exception as e:
                    print(f"Error reading {csv_file}: {e}")

    return data

def plot_metric(col_idx, title_label, ylabel, filename_suffix, data):
    plt.figure(figsize=(8, 5))

    # Use real MPI process counts on x-axis (not equally spaced)
    x = np.array(PROCESSES)

    for partition, color in PARTITIONS.items():
        y = [data[col_idx][partition].get(p, np.nan) for p in PROCESSES]
        plt.plot(x, y, marker="o", color=color, label=partition)

    plt.xlabel("Number of MPI Processes")
    plt.ylabel(ylabel)
    plt.title(f"{MATRIX_NAME} - {title_label}")
    plt.xticks(PROCESSES, PROCESSES)
    plt.legend()
    plt.tight_layout()

    # Horizontal reference lines (one for each y-tick)
    ax = plt.gca()
    for ytick in ax.get_yticks():
        plt.axhline(y=ytick, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

    out = os.path.join(PLOTS_DIR, f"{MATRIX_NAME}_{filename_suffix}(v2).png")
    plt.savefig(out)
    plt.close()

    print(f"Saved {out}")

def main():
    matrix_path = os.path.join(RESULTS_DIR, MATRIX_NAME)
    if not os.path.isdir(matrix_path):
        print(f"ERROR: matrix folder '{matrix_path}' not found.")
        return

    data = collect_data(matrix_path)

    for col_idx, title_label, ylabel, filename_suffix in METRICS:
        plot_metric(col_idx, title_label, ylabel, filename_suffix, data)

if __name__ == "__main__":
    main()
