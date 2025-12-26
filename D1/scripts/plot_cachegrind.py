import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

base_path = "../results/cachegrind_logs/"

# Patterns to capture references and misses
refs_patterns = {
    "I": re.compile(r"I\s+refs:\s+([\d,]+)"),
    "D": re.compile(r"D\s+refs:\s+([\d,]+)"),
    "LL": re.compile(r"LL\s+refs:\s+([\d,]+)")
}

miss_patterns = {
    "I1": re.compile(r"I1\s+misses:\s+([\d,]+)"),
    "LLi": re.compile(r"LLi\s+misses:\s+([\d,]+)"),
    "D1": re.compile(r"D1\s+misses:\s+([\d,]+)"),
    "LLd": re.compile(r"LLd\s+misses:\s+([\d,]+)"),
    "LL": re.compile(r"LL\s+misses:\s+([\d,]+)")
}

data = []

for matrix_folder in os.listdir(base_path):
    folder_path = os.path.join(base_path, matrix_folder)
    if not os.path.isdir(folder_path):
        continue

    err_file = os.path.join(folder_path, f"{matrix_folder}.err")
    if not os.path.isfile(err_file):
        continue

    with open(err_file, "r") as f:
        content = f.read()

    row = {"Matrix": matrix_folder}

    # Extract references
    refs = {}
    for key, pattern in refs_patterns.items():
        match = pattern.search(content)
        refs[key] = int(match.group(1).replace(",", "")) if match else 1  # avoid div by zero

    # Extract misses and calculate percentage
    for key, pattern in miss_patterns.items():
        match = pattern.search(content)
        if match:
            misses = int(match.group(1).replace(",", ""))
            # Determine which reference to use
            if key in ["I1", "LLi"]:
                ref = refs["I"]
            elif key in ["D1", "LLd"]:
                ref = refs["D"]
            elif key == "LL":
                ref = refs["LL"]
            else:
                ref = 1
            row[key] = (misses / ref) * 100
        else:
            row[key] = 0.0

    data.append(row)

df = pd.DataFrame(data)
df = df.set_index("Matrix")

# Plot grouped bar chart
cache_types = ["I1", "LLi", "D1", "LLd", "LL"]
x = np.arange(len(df.index))
width = 0.15

plt.figure(figsize=(12, 6))
for i, cache in enumerate(cache_types):
    plt.bar(x + i * width, df[cache], width=width, label=cache)

plt.xticks(x + width * 2, df.index, rotation=45)
plt.ylabel("Miss rate (%)")
plt.title("Cache Miss Rates per Matrix")
plt.legend()
plt.tight_layout()
plt.show()
