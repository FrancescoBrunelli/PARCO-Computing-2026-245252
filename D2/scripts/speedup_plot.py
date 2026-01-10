import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RESULTS_DIR = "results"
PLOTS_DIR = "plots"

os.makedirs(PLOTS_DIR, exist_ok=True)

# Partitionings
PARTITIONS = ["1D_Block", "1D_Cyclic", "2D_Block"]
COLORS = ["blue", "red", "green"]
BAR_WIDTH = 0.22

# Include NP = 1 (sequential baseline)
PROCESSES = [1, 2, 4, 8, 16, 32, 64, 128]

def read_p90_exec_time(csv_path):
    """Return p90 execution time (ms) from column 0."""
    df = pd.read_csv(csv_path, comment="#", header=None)
    return np.percentile(df[0].values, 90)

def plot_speedup(matrix):
    matrix_path = os.path.join(RESULTS_DIR, matrix)
    if not os.path.isdir(matrix_path):
        return

    # Sequential baseline
    seq_file = os.path.join(matrix_path, f"{matrix}_NP1.csv")
    if not os.path.isfile(seq_file):
        print(f"Skipping {matrix}: missing NP1 baseline")
        return

    seq_p90 = read_p90_exec_time(seq_file)

    speedup_data = {p: [] for p in PARTITIONS}

    for partition in PARTITIONS:
        part_path = os.path.join(matrix_path, partition)
        vals = []

        for np_ in PROCESSES:
            # Sequential case
            if np_ == 1:
                vals.append(1.0)
                continue

            csv_file = os.path.join(part_path, f"{matrix}_NP{np_}.csv")
            if not os.path.isfile(csv_file):
                vals.append(np.nan)
                continue

            par_p90 = read_p90_exec_time(csv_file)
            vals.append(seq_p90 / par_p90)

        speedup_data[partition] = vals

    # Numeric x positions (equally spaced)
    x = np.arange(len(PROCESSES), dtype=float)

    plt.figure(figsize=(12, 6))

    for i, partition in enumerate(PARTITIONS):
        offset = (i - 1) * BAR_WIDTH
        plt.bar(
            x + offset,
            speedup_data[partition],
            width=BAR_WIDTH,
            color=COLORS[i],
            label=partition,
            edgecolor="black"
        )

    plt.xticks(x, PROCESSES)
    plt.xlabel("Number of Processes")
    plt.ylabel("Speedup")
    plt.title(f"Strong Scaling — {matrix}")
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    plt.legend(title="Partitioning")
    plt.tight_layout()

    out_file = os.path.join(PLOTS_DIR, f"{matrix}_speedup.png")
    plt.savefig(out_file)
    plt.close()

    print(f"Saved: {out_file}")

def main():
    for matrix in os.listdir(RESULTS_DIR):
        plot_speedup(matrix)

if __name__ == "__main__":
    main()
