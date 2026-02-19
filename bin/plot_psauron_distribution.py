import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from statsmodels.nonparametric.kde import KDEUnivariate
import pandas as pd
import sys
import os

# =========================================================
# ARGUMENTS
# =========================================================
if len(sys.argv) < 2:
    print("Usage: python plot_psauron_distribution.py <combined_csv> [highlight_csv]")
    sys.exit(1)

csv_file = sys.argv[1]
highlight_csv = sys.argv[2] if len(sys.argv) >= 3 else None

# =========================================================
# LOAD CSV
# =========================================================
df = pd.read_csv(csv_file)
print("Loaded CSV")

species_to_plot = [
    "Mus_musculus",
    "Drosophila_melanogaster",
    "Arabidopsis_thaliana",
    "Saccharomyces_cerevisiae",
    "Rattus_norvegicus",
    "Homo_sapiens",
    "CaurisB8441",
    "Pan_troglodytes"
]

species_list = df["species"].unique()
print("Species found:", species_list)

# =========================================================
# LOAD HIGHLIGHT CSV (OPTIONAL)
# =========================================================
highlight_values = []
highlight_label = None

if highlight_csv and os.path.exists(highlight_csv):
    highlight_label = os.path.basename(highlight_csv).replace(".csv", "")
    print("Loading highlight CSV:", highlight_label)

    tmp_df = pd.read_csv(highlight_csv, skiprows=2)
    if "in-frame_score" in tmp_df.columns:
        highlight_values = tmp_df["in-frame_score"].dropna().to_numpy()

# =========================================================
# PREP DATA (FILTERED SPECIES ONLY)
# =========================================================
species_values = {}
for species, subdf in df.groupby("species"):
    if species not in species_to_plot:
        continue
    species_values[species] = subdf["in-frame_score"].dropna().to_numpy()

# =========================================================
# STYLE FUNCTIONS
# =========================================================
def style_background():
    return dict(linestyle="--", alpha=0.3, linewidth=1.5)

def style_highlight():
    return dict(linestyle="-", alpha=1.0, linewidth=2.5)

# =========================================================
# HISTOGRAM
# =========================================================
plt.figure()
for species, values in species_values.items():
    if len(values) == 0:
        continue
    plt.hist(values, bins=40, histtype="step", density=True,
             label=species, **style_background())

if len(highlight_values) > 0:
    plt.hist(highlight_values, bins=40, histtype="step", density=True,
             label=highlight_label, **style_highlight())

plt.title("Psauron In-frame Score Histogram")
plt.xlabel("In-frame Score")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig("psauron_hist.png", dpi=300)
plt.close()

# =========================================================
# STATS MODELS KDE
# =========================================================
plt.figure()
for species, values in species_values.items():
    if len(values) < 2:
        continue
    kde = KDEUnivariate(values)
    kde.fit(bw="scott")
    plt.plot(kde.support, kde.density,
             label=species, **style_background())

if len(highlight_values) > 1:
    kde = KDEUnivariate(highlight_values)
    kde.fit(bw="scott")
    plt.plot(kde.support, kde.density,
             label=highlight_label, **style_highlight())

plt.title("Psauron In-frame Score Density (StatsModels KDE)")
plt.xlabel("In-frame Score")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig("psauron_density_statsmodels.png", dpi=300)
plt.close()

# =========================================================
# SCIPY KDE
# =========================================================
plt.figure()
grid = np.linspace(0, 1, 400)

for species, values in species_values.items():
    if len(values) < 2:
        continue
    y = gaussian_kde(values)(grid)
    plt.plot(grid, y, label=species, **style_background())

if len(highlight_values) > 1:
    y = gaussian_kde(highlight_values)(grid)
    plt.plot(grid, y, label=highlight_label, **style_highlight())

plt.title("Psauron In-frame Score Density (SciPy KDE)")
plt.xlabel("In-frame Score")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig("psauron_density_scipy.png", dpi=300)
plt.close()

print("Done — selected species plotted with optional highlight overlay")

