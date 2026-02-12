import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from statsmodels.nonparametric.kde import KDEUnivariate
import sys
import glob
import pandas as pd
import os

if len(sys.argv) != 2:
    print("Usage: python plot_plddt.py <dataset_name>")
    sys.exit(1)

dataset_name = sys.argv[1]

# ---------- CSV FILE ----------
csv_file = os.path.join(os.path.dirname(__file__), "plddt_model_organisms.csv")

# Only plot these species from CSV
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

# ---------- LOAD CSV ----------
csv_species_values = {}
try:
    df_csv = pd.read_csv(csv_file)

    # Filter to required species only
    df_csv = df_csv[df_csv["Species"].isin(species_to_plot)]

    csv_species_values = (
        df_csv.groupby("Species")["Mean_pLDDT"]
        .apply(list)
        .to_dict()
    )

    print("Loaded CSV species:", list(csv_species_values.keys()))

except Exception as e:
    print("CSV not loaded:", e)

# ---------- LOAD PKL ----------
pkl_files = glob.glob(f"plddt_all_values_{dataset_name}*.pkl")

if len(pkl_files) == 0:
    print(f"ERROR: file not found: plddt_all_values_{dataset_name}*.pkl")
    sys.exit(1)

pkl_file = pkl_files[0]

with open(pkl_file, "rb") as f:
    species_plddt = pickle.load(f)

print("Loaded PKL species:", list(species_plddt.keys()))

def scipy_kde_curve(values, grid):
    kde = gaussian_kde(values)
    return kde.evaluate(grid)

plt.rcParams["figure.figsize"] = (10, 6)

# ---------- CONVERT PKL DICT ----------
species_values = {}
for species, protein_dict in species_plddt.items():
    species_values[species] = list(protein_dict.values())

# ---------- CREATE KDE GRID ----------
all_values_list = [vals for vals in species_values.values() if len(vals) > 0]

if len(csv_species_values) > 0:
    all_values_list.extend([vals for vals in csv_species_values.values() if len(vals) > 0])

all_values = np.concatenate(all_values_list)
grid = np.linspace(min(all_values), max(all_values), 400)

# =========================================================
# FIGURE 1 — Histogram
# =========================================================
plt.figure()

# PKL dataset
for species, values in species_values.items():
    plt.hist(values, bins=40, histtype="step", density=True, label=species)

# CSV overlay (transparent dashed)
for species, values in csv_species_values.items():
    plt.hist(values, bins=40, histtype="step", density=True,
             linestyle="--", alpha=0.3, label=f"{species} (CSV)")

plt.title(f"pLDDT Histogram ({dataset_name})")
plt.xlabel("Mean pLDDT")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig(f"plddt_hist_{dataset_name}.png", dpi=300)
plt.close()

# =========================================================
# FIGURE 2 — StatsModels KDE
# =========================================================
plt.figure()

# PKL dataset
for species, values in species_values.items():
    if len(values) < 2:
        continue
    kde = KDEUnivariate(values)
    kde.fit(bw="scott")
    plt.plot(kde.support, kde.density, label=species)

# CSV overlay
for species, values in csv_species_values.items():
    if len(values) < 2:
        continue
    kde = KDEUnivariate(values)
    kde.fit(bw="scott")
    plt.plot(kde.support, kde.density,
             linestyle="--", alpha=0.3, label=f"{species} (CSV)")

plt.title(f"pLDDT Density (StatsModels KDE) - {dataset_name}")
plt.xlabel("Mean pLDDT")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig(f"plddt_density_statsmodels_{dataset_name}.png", dpi=300)
plt.close()

# =========================================================
# FIGURE 3 — SciPy KDE
# =========================================================
plt.figure()

# PKL dataset
for species, values in species_values.items():
    if len(values) < 2:
        continue
    y = scipy_kde_curve(values, grid)
    plt.plot(grid, y, label=species)

# CSV overlay
for species, values in csv_species_values.items():
    if len(values) < 2:
        continue
    y = scipy_kde_curve(values, grid)
    plt.plot(grid, y,
             linestyle="--", alpha=0.3, label=f"{species} (CSV)")

plt.title(f"pLDDT Density (SciPy KDE) - {dataset_name}")
plt.xlabel("Mean pLDDT")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig(f"plddt_density_scipy_{dataset_name}.png", dpi=300)
plt.close()


