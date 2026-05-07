import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import sys
import glob
import pandas as pd
import os

from scipy.stats import skew
from sklearn.mixture import GaussianMixture
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from adjustText import adjust_text

if len(sys.argv) != 2:
    print("Usage: python plot_plddt.py <dataset_name>")
    sys.exit(1)

dataset_name = sys.argv[1]

# ---------- CSV FILE ----------
csv_file = os.path.join(os.path.dirname(__file__), "..", "reference", "plddt_model_organisms.csv")

# Cache file
cache_file = "cached_csv_gmm_metrics_model_organisms.csv"

# Species to exclude from CSV
species_to_exclude = ["Cauris6684", "Homo_sapiens_2k", "afdb_Homo_sapiens_2k"]

# ---------- LOAD CSV (all species except excluded, HS renamed) ----------
csv_species_values = {}
try:
    df_csv = pd.read_csv(csv_file)
    df_csv["Species"] = df_csv["Species"].replace({"HS": "Homo_sapiens"})
    df_csv = df_csv[~df_csv["Species"].isin(species_to_exclude)]
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

# =========================================================
# GUIDED GMM INITIALISATION from reference species
# =========================================================
reference_species = [
    "Drosophila_melanogaster",
    "Homo_sapiens",
    "Rattus_norvegicus",
    "Mus_musculus",
    "Pan_troglodytes",
    "Arabidopsis_thaliana"
]

print("Fitting reference species to compute initial GMM parameters...")

gmm_means_list = []
gmm_cov_list = []
gmm_weights_list = []

for sp in reference_species:
    if sp in species_values:
        values = np.array(species_values[sp])
    elif sp in csv_species_values:
        values = np.array(csv_species_values[sp])
    else:
        print(f"  Reference species not found, skipping: {sp}")
        continue

    if len(values) < 20:
        continue

    X = values.reshape(-1, 1)
    gmm = GaussianMixture(n_components=2, random_state=0).fit(X)

    means = gmm.means_.flatten()
    order = np.argsort(means)

    gmm_means_list.append(means[order])
    gmm_cov_list.append(gmm.covariances_.flatten()[order])
    gmm_weights_list.append(gmm.weights_[order])

if len(gmm_means_list) == 0:
    print("WARNING: No reference species found — falling back to default GMM init.")
    reference_means = None
    reference_weights = None
else:
    reference_means = np.mean(gmm_means_list, axis=0)
    reference_weights = np.mean(gmm_weights_list, axis=0)
    print(f"  Reference means:   {reference_means}")
    print(f"  Reference weights: {reference_weights}")

print("Reference GMM parameters computed.")

def fit_gmm(values):
    """Fit 2-component GMM with guided init if available, else default."""
    X = values.reshape(-1, 1)
    if reference_means is not None:
        gmm = GaussianMixture(
            n_components=2,
            means_init=reference_means.reshape(-1, 1),
            weights_init=reference_weights,
            random_state=0
        ).fit(X)
    else:
        gmm = GaussianMixture(n_components=2, random_state=0).fit(X)
    return gmm

# ---------- CREATE KDE GRID ----------
all_values_list = [vals for vals in species_values.values() if len(vals) > 0]
if len(csv_species_values) > 0:
    all_values_list.extend([vals for vals in csv_species_values.values() if len(vals) > 0])

all_values = np.concatenate(all_values_list)
grid = np.linspace(min(all_values), max(all_values), 400)

# =========================================================
# FIGURE — SciPy KDE density plot
# =========================================================
plt.figure()

for species, values in species_values.items():
    if len(values) < 2:
        continue
    y = scipy_kde_curve(values, grid)
    plt.plot(grid, y, label=species)

for species, values in csv_species_values.items():
    if len(values) < 2:
        continue
    y = scipy_kde_curve(values, grid)
    plt.plot(grid, y, linestyle="--", alpha=0.3, label=f"{species} (CSV)")

plt.title(f"pLDDT Density (SciPy KDE) - {dataset_name}")
plt.xlabel("Mean pLDDT")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig(f"plddt_density_scipy_{dataset_name}.png", dpi=300)
plt.close()

# =========================================================
# GMM SCATTER
# =========================================================
metrics = []

# ---------------------------
# LOAD OR COMPUTE CSV METRICS
# ---------------------------
if os.path.exists(cache_file):
    print("Loading cached CSV GMM metrics...")
    csv_metrics_df = pd.read_csv(cache_file)
    # Apply rename and exclusions to cache in case it was built with old names
    csv_metrics_df["Species"] = csv_metrics_df["Species"].replace({"HS": "Homo_sapiens"})
    csv_metrics_df = csv_metrics_df[~csv_metrics_df["Species"].isin(species_to_exclude)]
else:
    print("Computing CSV GMM metrics (first run)...")
    csv_metrics = []

    for sp, values in csv_species_values.items():
        values = np.array(values)
        if len(values) < 20:
            continue
        gmm = fit_gmm(values)
        means = gmm.means_.flatten()
        order = np.argsort(means)
        upper_prop = gmm.weights_[order][1]
        csv_metrics.append({
            "Species": sp,
            "Skew": skew(values),
            "GMM_upper_prop": upper_prop,
            "Prop_ge_70": np.mean(values >= 70)
        })

    csv_metrics_df = pd.DataFrame(csv_metrics)
    csv_metrics_df.to_csv(cache_file, index=False)
    print("CSV metrics cached.")

for _, row in csv_metrics_df.iterrows():
    metrics.append({
        "Species": row["Species"],
        "Skew": row["Skew"],
        "GMM_upper_prop": row["GMM_upper_prop"],
        "Prop_ge_70": row["Prop_ge_70"],
        "Is_PKL": False
    })

# ---------------------------
# PKL METRICS (always fresh)
# ---------------------------
for sp, values in species_values.items():
    values = np.array(values)
    if len(values) < 20:
        continue
    gmm = fit_gmm(values)
    means = gmm.means_.flatten()
    order = np.argsort(means)
    upper_prop = gmm.weights_[order][1]
    metrics.append({
        "Species": sp,
        "Skew": skew(values),
        "GMM_upper_prop": upper_prop,
        "Prop_ge_70": np.mean(values >= 70),
        "Is_PKL": True
    })

df_metrics = pd.DataFrame(metrics)

# ---------------------------
# SCATTER — all points coloured by skew, PKL query highlighted
# ---------------------------
norm = mcolors.Normalize(
    vmin=df_metrics["Skew"].min(),
    vmax=df_metrics["Skew"].max()
)
cmap = cm.viridis

fig, ax = plt.subplots(figsize=(9, 7))

# Plot all points together coloured by skew
sc = ax.scatter(
    df_metrics["GMM_upper_prop"],
    df_metrics["Prop_ge_70"],
    c=df_metrics["Skew"],
    cmap=cmap,
    norm=norm,
    alpha=0.8,
    zorder=2
)

fig.colorbar(sc, ax=ax, label="Skew")

# Expand axis limits by 10% on each side to prevent clipping
x_vals = df_metrics["GMM_upper_prop"]
y_vals = df_metrics["Prop_ge_70"]
x_margin = (x_vals.max() - x_vals.min()) * 0.10
y_margin = (y_vals.max() - y_vals.min()) * 0.10
ax.set_xlim(x_vals.min() - x_margin, x_vals.max() + x_margin)
ax.set_ylim(y_vals.min() - y_margin, y_vals.max() + y_margin)

# Build text objects for adjustText — CSV species
texts = []
for _, row in df_metrics[~df_metrics["Is_PKL"]].iterrows():
    t = ax.text(
        row["GMM_upper_prop"],
        row["Prop_ge_70"],
        row["Species"],
        fontsize=8
    )
    texts.append(t)

# Overlay PKL query species with red circle — no legend entry
pkl_points = df_metrics[df_metrics["Is_PKL"]]

if len(pkl_points) > 0:
    ax.scatter(
        pkl_points["GMM_upper_prop"],
        pkl_points["Prop_ge_70"],
        facecolors="none",
        edgecolors="red",
        s=200,
        linewidths=2,
        zorder=3
    )
    for _, row in pkl_points.iterrows():
        t = ax.text(
            row["GMM_upper_prop"],
            row["Prop_ge_70"],
            row["Species"],
            fontsize=9,
            fontweight="bold",
            color="red"
        )
        texts.append(t)

# Auto-adjust all labels to avoid overlap
adjust_text(
    texts,
    ax=ax,
    expand=(1.3, 1.5),
    arrowprops=dict(arrowstyle="-", color="grey", lw=0.5)
)

ax.set_xlabel("GMM Upper Component Proportion")
ax.set_ylabel("Proportion ≥70 pLDDT")
ax.set_title(f"GMM Upper Peak vs ≥70 Proportion ({dataset_name})")
plt.tight_layout()
plt.savefig(f"plddt_gmm_scatter_{dataset_name}.png", dpi=300)
plt.close()

print("GMM scatter diagnostics generated (CSV cached, PKL highlighted).")
