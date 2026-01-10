import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Base results folder
RESULTS_DIR = "results"
PLOTS_DIR = "plots"

# Ensure plots folder exists
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

# Processes (x-axis)
PROCESSES = [2, 4, 8, 16, 32, 64, 128]

def extract_csv_p90(filepath, col_idx, comm_time=False):
    df = pd.read_csv(filepath, comment="#", header=None)
    values = df[col_idx].values
    p90 = np.percentile(values, 90)
    if comm_time:
        p90 /= 1000
    return p90

def collect_data(matrix_path):
    data = {col_idx: {p: {} for p in PROCESSES} for col_idx, _, _, _ in METRICS}

    for partition, color in PARTITIONS.items():
        part_path = os.path.join(matrix_path, partition)
        if not os.path.isdir(part_path):
            continue

        for p in PROCESSES:
            csv_file = os.path.join(
                part_path, f"{os.path.basename(matrix_path)}_NP{p}.csv"
            )
            if not os.path.isfile(csv_file):
                print(f"Warning: missing {csv_file}")
                continue

            for col_idx, _, _, _ in METRICS:
                comm_time = (col_idx == 5)
                try:
                    val = extract_csv_p90(csv_file, col_idx, comm_time)
                    data[col_idx][p][partition] = val
                except Exception as e:
                    print(f"Error reading {csv_file}: {e}")

    return data

def plot_metric(matrix_name, col_idx, title_label, ylabel, filename_suffix, data):
    plt.figure(figsize=(8, 5))

    # Use indices for equally spaced plotting
    x_indices = np.arange(len(PROCESSES))

    for partition, color in PARTITIONS.items():
        y = [data[col_idx][p].get(partition, np.nan) for p in PROCESSES]
        plt.plot(x_indices, y, color=color, marker="o", label=partition)

    # X-axis
    plt.xticks(x_indices, PROCESSES)
    plt.xlabel("Number of Processes")

    # Labels and title
    plt.ylabel(ylabel)
    plt.title(f"{matrix_name} - {title_label}")
    plt.legend()
    plt.tight_layout()

    # Horizontal reference lines (one per y-tick)
    ax = plt.gca()
    for y in ax.get_yticks():
        plt.axhline(y=y, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

    # Save plot
    plot_filename = f"{matrix_name}_{filename_suffix}.png"
    plt.savefig(os.path.join(PLOTS_DIR, plot_filename))
    plt.close()

def main():
    for matrix in os.listdir(RESULTS_DIR):
        matrix_path = os.path.join(RESULTS_DIR, matrix)
        if not os.path.isdir(matrix_path):
            continue

        print(f"Processing matrix: {matrix}")
        data = collect_data(matrix_path)

        for col_idx, title_label, ylabel, filename_suffix in METRICS:
            plot_metric(matrix, col_idx, title_label, ylabel, filename_suffix, data)

if __name__ == "__main__":
    main()
