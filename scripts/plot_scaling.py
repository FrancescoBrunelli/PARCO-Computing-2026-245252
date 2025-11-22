import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

if len(sys.argv) < 2:
    print("Usage: python3 plot_speedup_hist.py <matrix_name>")
    sys.exit(1)

MATRIX = sys.argv[1]
CSV_FILE = os.path.join("results", "runtime_summary.csv")

MEMORY = "64GB"
SCHEDULINGS = ["static", "dynamic", "guided"]
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c"]
BAR_WIDTH = 0.22

df = pd.read_csv(CSV_FILE)

df = df[(df["matrix"] == MATRIX) & (df["memory"] == MEMORY)]
if df.empty:
    print(f"ERROR: no data found for matrix '{MATRIX}' with memory {MEMORY}.")
    sys.exit(1)

# Sort thread levels
thread_levels = sorted(df["threads"].unique())
x = np.arange(len(thread_levels))  # base positions

plt.figure(figsize=(12, 6))

for i, sched in enumerate(SCHEDULINGS):
    df_s = df[df["scheduling"] == sched].sort_values("threads")

    if df_s.empty:
        continue

    seq_time = df_s[df_s["threads"] == 1]["p90_runtime_ms"].values[0]
    speedups = seq_time / df_s["p90_runtime_ms"]

    offset = (i - 1) * BAR_WIDTH   # centers the 3 bars around x

    plt.bar(
        x + offset, speedups,
        width=BAR_WIDTH,
        color=COLORS[i],
        label=sched,
        edgecolor='black'
    )

plt.xticks(x, thread_levels)
plt.xlabel("Number of Threads")
plt.ylabel("Speedup")
plt.title(f"Strong Scaling — {MATRIX}")

plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.legend(title="Scheduling")
plt.tight_layout()

plt.show()
