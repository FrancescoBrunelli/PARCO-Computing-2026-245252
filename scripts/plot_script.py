import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
ram_sizes = ["32GB", "64GB", "128GB"]
data = {}

for ram in ram_sizes:
    x_list, y_list = [], []
    directory = f"../results/epb1/{ram}/runtime/"
    for file in os.scandir(directory):
        with open(os.path.join(directory, file.name)) as f:
            threads = int(f.readline().strip())
            df = pd.read_csv(f, comment="#", skiprows=1)
            p90 = np.percentile(df.iloc[:,0].values, 90)
            x_list.append(threads)
            y_list.append(p90)
            print(f"{ram} - File {file.name}: p90 = {p90}")
    # Sort
    x_list, y_list = zip(*sorted(zip(x_list, y_list)))
    data[ram] = (x_list, y_list)

# Plot
threads_ticks = data["32GB"][0]  # x-values for 32GB, same for all
plt.xscale("log")
plt.xticks(threads_ticks, threads_ticks)
#plt.ylim(0.42, 0.48)

plt.grid(True, which="both", axis='y', ls='--')

colors = {"32GB":"blue", "64GB":"green", "128GB":"red"}
for ram in ram_sizes:
    plt.plot(*data[ram], color=colors[ram], label=ram)

plt.title("epb1")
plt.xlabel("No. of Threads")
plt.ylabel("Execution Time")
plt.legend(title = "Memory")
plt.show()
