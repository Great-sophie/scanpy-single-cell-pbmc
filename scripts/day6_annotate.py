from pathlib import Path
import scanpy as sc
import pandas as pd

OUT_DIR = Path("outputs")
FIG_DIR = Path("figures")

sc.settings.figdir = str(FIG_DIR)

print("=== Load processed PBMC3k ===")

adata = sc.read_h5ad(
    OUT_DIR / "pbmc3k_processed.h5ad"
)

print(adata)

# --------------------------------------------------
# 1. Cluster -> cell type mapping
# --------------------------------------------------

cluster_to_celltype = {
    "0": "T cells",
    "1": "CD14-like monocytes",
    "2": "NK / cytotoxic lymphocytes",
    "3": "B cells",
    "4": "Dendritic cells",
    "5": "Platelets",
    "6": "Cycling cells",
}

adata.obs["cell_type"] = (
    adata.obs["leiden"]
    .astype(str)
    .map(cluster_to_celltype)
)

print("\n=== Cell type counts ===")
print(
    adata.obs["cell_type"]
    .value_counts()
)

# --------------------------------------------------
# 2. Cell-type UMAP
# --------------------------------------------------

sc.pl.umap(
    adata,
    color="cell_type",
    legend_loc="right margin",
    title="PBMC3k cell-type annotation",
    save="_08_cell_types.png",
)

# --------------------------------------------------
# 3. Canonical marker DotPlot
# --------------------------------------------------

marker_dict = {
    "T cells": [
        "CD3D",
        "CD3E",
        "IL7R",
        "CCR7",
    ],

    "Cytotoxic / NK": [
        "NKG7",
        "GNLY",
        "PRF1",
        "CCL5",
    ],

    "B cells": [
        "MS4A1",
        "CD79A",
        "CD74",
    ],

    "Monocytes": [
        "LYZ",
        "S100A8",
        "S100A9",
        "LST1",
    ],

    "Dendritic": [
        "FCER1A",
        "CST3",
        "HLA-DRA",
    ],

    "Platelets": [
        "PPBP",
        "PF4",
    ],

    "Cycling": [
        "TYMS",
        "PCNA",
    ],
}

available_marker_dict = {}

for group, genes in marker_dict.items():
    available = [
        g for g in genes
        if g in adata.raw.var_names
    ]

    if available:
        available_marker_dict[group] = available

sc.pl.dotplot(
    adata,
    available_marker_dict,
    groupby="cell_type",
    use_raw=True,
    standard_scale="var",
    save="_09_celltype_marker_dotplot.png",
)

# --------------------------------------------------
# 4. Cell-type composition
# --------------------------------------------------

summary = (
    adata.obs["cell_type"]
    .value_counts()
    .rename_axis("cell_type")
    .reset_index(name="n_cells")
)

summary["fraction"] = (
    summary["n_cells"]
    / summary["n_cells"].sum()
)

summary.to_csv(
    OUT_DIR / "cell_type_summary.csv",
    index=False,
)

print("\n=== Cell-type composition ===")
print(summary)

# --------------------------------------------------
# 5. Save annotated AnnData
# --------------------------------------------------

adata.write(
    OUT_DIR / "pbmc3k_annotated.h5ad"
)

adata.obs.to_csv(
    OUT_DIR / "cell_metadata_annotated.csv"
)

print("\nSaved:")
print(
    OUT_DIR / "pbmc3k_annotated.h5ad"
)

print("\nDay 6 annotation complete.")
