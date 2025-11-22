import os
import sys
import numpy as np
import matplotlib.pyplot as plt

RESULTS_DIR = "results"
RAM_SIZES = ["32GB", "64GB", "128GB"]
SCHEDULINGS = ["static", "dynamic", "guided"]

COLORS = {"32GB": "blue", "64GB": "green", "128GB": "red"}

def extract_times(filepath):
    times = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Skip thread count (first numeric line)
            if line.isdigit():
                continue
            try:
                times.append(float(line))
            except ValueError:
                pass
    return times

def collect_data(matrix_name, scheduling):
    data = {}

    for ram in RAM_SIZES:
        x_list, y_list = [], []
        directory = os.path.join(RESULTS_DIR, matrix_name, scheduling, ram, "runtime")

        if not os.path.isdir(directory):
            print(f"Missing directory: {directory}")
            continue

        for file in sorted(os.listdir(directory)):
            if not file.endswith(".txt"):
                continue

            filepath = os.path.join(directory, file)
            times = extract_times(filepath)
            if not times:
                continue

            threads = None
            for token in file.replace(".txt", "").split("_"):
                if token.startswith("NT"):
                    threads = int(token[2:])
                    break
            if threads is None:
                # fallback: read from file content
                with open(filepath) as f:
                    for line in f:
                        if line.strip().isdigit():
                            threads = int(line.strip())
                            break

            if threads is None:
                continue

            p90 = np.percentile(times, 90)
            x_list.append(threads)
            y_list.append(p90)
            print(f"{scheduling.upper()} | {ram} | {file}: p90 = {p90:.6f}")

        if x_list:
            x_list, y_list = zip(*sorted(zip(x_list, y_list)))
            data[ram] = (x_list, y_list)

    return data

def plot_data(matrix_name, scheduling, data):
    #Create a plot for a specific scheduling type.
    if not data:
        print(f"No data available for scheduling '{scheduling}'")
        return

    plt.figure(figsize=(8, 5))
    plt.xscale("log")
    plt.grid(True, which="both", axis="y", ls="--")

    threads_ticks = next(iter(data.values()))[0]
    plt.xticks(threads_ticks, threads_ticks)

    for ram, (x, y) in data.items():
        plt.plot(x, y, color=COLORS.get(ram, "black"), marker="o", label=ram)

    plt.title(f"{matrix_name} ({scheduling})")
    plt.xlabel("No. of Threads")
    plt.ylabel("Execution Time [ms]")
    plt.legend(title="Memory")
    plt.tight_layout()
    plt.show()

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_runtime_all_sched.py <matrix_name>")
        sys.exit(1)

    matrix_name = sys.argv[1]

    for sched in SCHEDULINGS:
        print(f"\n=== Processing scheduling type: {sched.upper()} ===")
        data = collect_data(matrix_name, sched)
        plot_data(matrix_name, sched, data)

if __name__ == "__main__":
    main()
