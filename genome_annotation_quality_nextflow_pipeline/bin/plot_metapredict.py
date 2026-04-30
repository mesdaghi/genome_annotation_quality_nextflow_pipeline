import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from statsmodels.nonparametric.kde import KDEUnivariate
import pandas as pd
import sys
import os

# =========================================================
# CONFIG
# =========================================================
REFERENCE_PKL = os.path.join(os.path.dirname(__file__), "combined_proteome_disorder.pkl")

# =========================================================
# ARGUMENTS
# =========================================================
if len(sys.argv) < 2:
    print("Usage: python plot_metapredict.py <metapredict_csv>")
    sys.exit(1)

highlight_csv = sys.argv[1]
dataset_label = os.path.basename(highlight_csv).replace(".csv", "")

print("Dataset:", dataset_label)

# =========================================================
# LOAD REFERENCE PKL
# =========================================================
if not os.path.exists(REFERENCE_PKL):
    raise FileNotFoundError(f"Reference PKL not found: {REFERENCE_PKL}")

df = pd.read_pickle(REFERENCE_PKL)
print("Loaded reference PKL")

species_list = df["species"].unique()
print("Species found:", species_list)

# =========================================================
# LOAD HIGHLIGHT CSV
# =========================================================
highlight_mean_values = []

if os.path.exists(highlight_csv):
    print("Loading highlight CSV:", highlight_csv)

    with open(highlight_csv) as f:
        for line in f:
            parts = line.strip().split(",")
            try:
                scores = [float(x) for x in parts[2:] if x]
            except:
                continue
            if scores:
                highlight_mean_values.append(np.mean(scores))
else:
    print("WARNING: highlight CSV not found")

# =========================================================
# PREP DATA
# =========================================================
species_mean_values = {}
for species, subdf in df.groupby("species"):
    species_mean_values[species] = subdf["mean_disorder"].to_numpy()

# =========================================================
# STYLE FUNCTIONS
# =========================================================
def style_background():
    return dict(linestyle="--", alpha=0.3, linewidth=1.5)

def style_highlight():
    return dict(linestyle="-", alpha=1.0, linewidth=2.5)

# =========================================================
# SCIPY KDE
# =========================================================
plt.figure()
grid = np.linspace(0, 1, 400)

for species, values in species_mean_values.items():
    if len(values) < 2:
        continue
    y = gaussian_kde(values)(grid)
    plt.plot(grid, y, label=species, **style_background())

if highlight_mean_values and len(highlight_mean_values) > 1:
    y = gaussian_kde(highlight_mean_values)(grid)
    plt.plot(grid, y, label=dataset_label, **style_highlight())

plt.title("Mean Disorder Density (SciPy KDE)")
plt.xlabel("Mean Disorder")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig(f"{dataset_label}_mean_disorder_density_scipy.png", dpi=300)
plt.close()

print("Done — overlay plots generated")

