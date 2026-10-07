from pathlib import Path
import scanpy as sc

ROOT = Path(".")
FIG_DIR = ROOT / "figures"
OUT_DIR = ROOT / "outputs"

FIG_DIR.mkdir(exist_ok=True)
OUT_DIR.mkdir(exist_ok=True)

sc.settings.verbosity = 3
sc.settings.figdir = str(FIG_DIR)

print("\n=== 1. Load PBMC3k ===")
adata = sc.datasets.pbmc3k()

print(adata)
print("Cells:", adata.n_obs)
print("Genes:", adata.n_vars)

# --------------------------------------------------
# 2. QC
# --------------------------------------------------

print("\n=== 2. QC metrics ===")

adata.var["mt"] = (
    adata.var_names
    .str.upper()
    .str.startswith("MT-")
)

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

print(
    adata.obs[
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
        ]
    ].describe()
)

sc.pl.violin(
    adata,
    [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
    ],
    jitter=0.4,
    multi_panel=True,
    save="_01_qc_violin.png",
)

# --------------------------------------------------
# 3. Filtering
# --------------------------------------------------

print("\n=== 3. Filtering ===")

before_cells = adata.n_obs
before_genes = adata.n_vars

sc.pp.filter_cells(
    adata,
    min_genes=200,
)

sc.pp.filter_genes(
    adata,
    min_cells=3,
)

adata = adata[
    adata.obs["pct_counts_mt"] < 5
].copy()

print(
    f"Cells: {before_cells} -> {adata.n_obs}"
)

print(
    f"Genes: {before_genes} -> {adata.n_vars}"
)

# --------------------------------------------------
# 4. Normalization
# --------------------------------------------------

print("\n=== 4. Normalize ===")

sc.pp.normalize_total(
    adata,
    target_sum=1e4,
)

sc.pp.log1p(adata)

adata.raw = adata

# --------------------------------------------------
# 5. HVG
# --------------------------------------------------

print("\n=== 5. Highly variable genes ===")

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
)

print(
    "HVGs:",
    int(adata.var["highly_variable"].sum())
)

sc.pl.highly_variable_genes(
    adata,
    save="_02_hvg.png",
)

adata = adata[
    :,
    adata.var["highly_variable"],
].copy()

# --------------------------------------------------
# 6. Scale
# --------------------------------------------------

print("\n=== 6. Scale ===")

sc.pp.scale(
    adata,
    max_value=10,
)

# --------------------------------------------------
# 7. PCA
# --------------------------------------------------

print("\n=== 7. PCA ===")

sc.tl.pca(
    adata,
    svd_solver="arpack",
)

sc.pl.pca_variance_ratio(
    adata,
    log=True,
    save="_03_pca_variance.png",
)

# --------------------------------------------------
# 8. Neighbors
# --------------------------------------------------

print("\n=== 8. Neighbor graph ===")

sc.pp.neighbors(
    adata,
    n_neighbors=10,
    n_pcs=40,
)

# --------------------------------------------------
# 9. UMAP
# --------------------------------------------------

print("\n=== 9. UMAP ===")

sc.tl.umap(adata)

sc.pl.umap(
    adata,
    save="_04_umap_unclustered.png",
)

# --------------------------------------------------
# 10. Leiden
# --------------------------------------------------

print("\n=== 10. Leiden clustering ===")

sc.tl.leiden(
    adata,
    resolution=0.5,
    key_added="leiden",
)

print(
    adata.obs["leiden"]
    .value_counts()
    .sort_index()
)

sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    save="_05_umap_leiden.png",
)

# --------------------------------------------------
# 11. Marker genes
# --------------------------------------------------

print("\n=== 11. Marker genes ===")

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    use_raw=True,
)

sc.pl.rank_genes_groups(
    adata,
    n_genes=15,
    sharey=False,
    save="_06_cluster_markers.png",
)

marker_df = sc.get.rank_genes_groups_df(
    adata,
    group=None,
)

marker_df.to_csv(
    OUT_DIR / "cluster_marker_genes.csv",
    index=False,
)

print(marker_df.head(30))

# --------------------------------------------------
# 12. Known immune markers
# --------------------------------------------------

print("\n=== 12. Known immune markers ===")

marker_genes = [
    "IL7R",
    "CCR7",
    "LTB",
    "CD3D",
    "CD8A",
    "NKG7",
    "GNLY",
    "MS4A1",
    "CD79A",
    "FCGR3A",
    "LYZ",
    "S100A8",
    "S100A9",
    "FCER1A",
    "CST3",
    "PPBP",
]

existing_markers = [
    g
    for g in marker_genes
    if g in adata.raw.var_names
]

sc.pl.umap(
    adata,
    color=existing_markers,
    use_raw=True,
    ncols=4,
    save="_07_known_markers.png",
)

# --------------------------------------------------
# 13. Save
# --------------------------------------------------

print("\n=== 13. Save outputs ===")

adata.obs.to_csv(
    OUT_DIR / "cell_metadata.csv"
)

cluster_summary = (
    adata.obs["leiden"]
    .value_counts()
    .sort_index()
    .rename_axis("cluster")
    .reset_index(name="n_cells")
)

cluster_summary["fraction"] = (
    cluster_summary["n_cells"]
    / cluster_summary["n_cells"].sum()
)

cluster_summary.to_csv(
    OUT_DIR / "cluster_summary.csv",
    index=False,
)

adata.write(
    OUT_DIR / "pbmc3k_processed.h5ad"
)

print(cluster_summary)

print("\nSaved:")
print(OUT_DIR / "pbmc3k_processed.h5ad")

print("\nDay 6 pipeline complete.")
