import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Configuration
CSV_FILE = "results/perf_summary.csv"
THREADS_HARDCODED = 64   # Change as needed
MEMORY_HARDCODED = "128GB"  # Change as needed

SCHEDULINGS = ["static", "dynamic", "guided"]
CACHE_METRICS = ["cache_misses_percent", "l1_misses_percent", "llc_misses_percent"]
MATRIX_ORDER = ["fd12", "sinc12", "c-46", "epb1", "human_gene2"]
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c"]
MIN_VISIBLE = 0.5
BAR_WIDTH = 0.2
CLUSTER_GAP = 0.6

# Load CSV
df = pd.read_csv(CSV_FILE)
df = df[(df["threads"] == THREADS_HARDCODED) & (df["memory"] == MEMORY_HARDCODED)]

fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)

for idx, sched in enumerate(SCHEDULINGS):
    df_sched = df[df["scheduling"] == sched].copy().set_index("matrix").reindex(MATRIX_ORDER)
    n_matrices = len(MATRIX_ORDER)
    n_metrics = len(CACHE_METRICS)
    cluster_width = n_metrics * BAR_WIDTH
    x = np.arange(n_matrices) * (cluster_width + CLUSTER_GAP)
    ax = axes[idx]

    for i, metric in enumerate(CACHE_METRICS):
        values = df_sched[metric].values.copy()
        mask_small = values < MIN_VISIBLE
        values_plot = np.where(mask_small, MIN_VISIBLE, values)
        offset = (i - (n_metrics - 1) / 2) * BAR_WIDTH
        bars = ax.bar(x + offset, values_plot, BAR_WIDTH, label=metric, color=COLORS[i])

        for j, bar in enumerate(bars):
            if mask_small[j]:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        MIN_VISIBLE + 0.5,
                        f"{values[j]:.2f}%",
                        ha='center', va='bottom', fontsize=3)
            else:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        values_plot[j] + 0.5,
                        f"{values[j]:.2f}%",
                        ha='center', va='bottom', fontsize=3)

    ax.set_xticks(x)
    ax.set_xticklabels(MATRIX_ORDER)
    ax.set_title(f"Scheduling: {sched}")
    ax.set_xlabel("Matrix")
    if idx == 0:
        ax.set_ylabel("Cache Miss Percentage [%]")
    ax.set_ylim(0, 100)
    ax.grid(axis="y", linestyle="--", alpha=0.7)

axes[-1].legend(CACHE_METRICS, loc='upper right')

plt.suptitle(f"Cache Misses per Matrix - Threads: {THREADS_HARDCODED}, Memory: {MEMORY_HARDCODED}")
plt.tight_layout(rect=[0, 0, 1, 0.95])
title = f"cache_misses_{MEMORY_HARDCODED}_{THREADS_HARDCODED}T"
downloads = os.path.join(os.path.expanduser("~"), "Downloads")
filepath = os.path.join(downloads, title)
plt.savefig(filepath, dpi=600)
plt.show()

# SINGLE PLOT VERSION
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

CSV_FILE = "results/perf_summary.csv"
THREADS_HARDCODED = 1   # Change as needed
MEMORY_HARDCODED = "128GB"  # Change as needed

SCHEDULINGS = ["static", "dynamic", "guided"]
CACHE_METRICS = ["cache_misses_percent", "l1_misses_percent", "llc_misses_percent"]
MATRIX_ORDER = ["fd12", "sinc12", "c-46", "epb1", "human_gene2"]
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c"]
MIN_VISIBLE = 0.5
BAR_WIDTH = 0.2
CLUSTER_GAP = 0.6

# Load CSV
df = pd.read_csv(CSV_FILE)
df = df[(df["threads"] == THREADS_HARDCODED) & (df["memory"] == MEMORY_HARDCODED)]

def plot_cache_misses(df_plot, title_suffix=""):
    n_matrices = len(MATRIX_ORDER)
    n_metrics = len(CACHE_METRICS)
    cluster_width = n_metrics * BAR_WIDTH
    x = np.arange(n_matrices) * (cluster_width + CLUSTER_GAP)

    fig, ax = plt.subplots(figsize=(12, 6))

    for i, metric in enumerate(CACHE_METRICS):
        values = df_plot[metric].values.copy()
        mask_small = values < MIN_VISIBLE
        values_plot = np.where(mask_small, MIN_VISIBLE, values)
        offset = (i - (n_metrics - 1) / 2) * BAR_WIDTH
        bars = ax.bar(x + offset, values_plot, BAR_WIDTH, label=metric, color=COLORS[i])

        for j, bar in enumerate(bars):
            if mask_small[j]:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        MIN_VISIBLE + 0.5,
                        f"<0.5% ({values[j]:.2f}%)",
                        ha='center', va='bottom', fontsize=7)
            else:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        values_plot[j] + 0.5,
                        f"{values[j]:.2f}%",
                        ha='center', va='bottom', fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(MATRIX_ORDER)
    ax.set_ylabel("Cache Miss Percentage [%]")
    ax.set_xlabel("Matrix")
    ax.set_title(f"Cache Misses per Matrix {title_suffix}")
    ax.set_ylim(0, 100)
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.7)

    plt.tight_layout()
    return fig, ax

for sched in SCHEDULINGS:
    df_sched = df[df["scheduling"] == sched].copy().set_index("matrix").reindex(MATRIX_ORDER)
    plot_cache_misses(df_sched, title_suffix=f"- Sequential, Memory: {MEMORY_HARDCODED}")
    plt.show()
'''