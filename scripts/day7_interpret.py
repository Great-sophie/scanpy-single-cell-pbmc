from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


OUT_DIR = Path("outputs")
FIG_DIR = Path("figures")

OUT_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

sc.settings.figdir = str(FIG_DIR)
sc.settings.verbosity = 2


# ============================================================
# 1. Load Day 6 annotated object
# ============================================================

print("\n=== 1. Load annotated PBMC3k ===")

adata = sc.read_h5ad(
    OUT_DIR / "pbmc3k_annotated.h5ad"
)

print(adata)
print("\nCurrent cell types:")
print(adata.obs["cell_type"].value_counts())


# ============================================================
# 2. Build lymphoid subtype scores using full normalized genes
# ============================================================

print("\n=== 2. Lymphoid signature scoring ===")

# .raw contains the full normalized/log1p gene matrix
raw = adata.raw.to_adata()
raw.obs = adata.obs.copy()

signatures = {
    "CD4_T": [
        "IL7R",
        "CCR7",
        "LTB",
        "MAL",
        "LTST1",
    ],
    "CD8_T": [
        "CD3D",
        "CD3E",
        "CD8A",
        "CCL5",
        "GZMK",
    ],
    "NK": [
        "NKG7",
        "GNLY",
        "PRF1",
        "FCGR3A",
        "TYROBP",
    ],
}

# Remove genes not present
for name, genes in signatures.items():

    genes_present = [
        g for g in genes
        if g in raw.var_names
    ]

    print(
        f"{name}: {genes_present}"
    )

    sc.tl.score_genes(
        raw,
        gene_list=genes_present,
        score_name=f"{name}_score",
        random_state=42,
    )

score_columns = [
    "CD4_T_score",
    "CD8_T_score",
    "NK_score",
]

for col in score_columns:
    adata.obs[col] = raw.obs[col].values


# ============================================================
# 3. Re-cluster lymphoid compartment
# ============================================================

print("\n=== 3. Lymphoid subclustering ===")

lymphoid_mask = adata.obs["cell_type"].isin(
    [
        "T cells",
        "NK / cytotoxic lymphocytes",
    ]
)

lym = adata[lymphoid_mask].copy()

print(
    "Lymphoid cells:",
    lym.n_obs
)

sc.tl.pca(
    lym,
    svd_solver="arpack",
)

n_pcs = min(
    30,
    lym.obsm["X_pca"].shape[1]
)

sc.pp.neighbors(
    lym,
    n_neighbors=10,
    n_pcs=n_pcs,
)

sc.tl.umap(
    lym,
    random_state=42,
)

sc.tl.leiden(
    lym,
    resolution=0.5,
    key_added="lymphoid_leiden",
    random_state=42,
)

print("\nLymphoid clusters:")
print(
    lym.obs["lymphoid_leiden"]
    .value_counts()
    .sort_index()
)


# ============================================================
# 4. Assign CD4-like / CD8-like / NK-like labels
# ============================================================

print("\n=== 4. Refine lymphoid annotation ===")

cluster_scores = (
    lym.obs
    .groupby("lymphoid_leiden", observed=True)
    [score_columns]
    .mean()
)

print("\nMean signature scores per lymphoid cluster:")
print(cluster_scores)

label_map = {
    "CD4_T_score": "CD4 T-like",
    "CD8_T_score": "CD8 T-like",
    "NK_score": "NK-like",
}

cluster_to_subtype = {}

for cluster, row in cluster_scores.iterrows():

    best_score = row.idxmax()

    cluster_to_subtype[str(cluster)] = (
        label_map[best_score]
    )

print("\nCluster → subtype:")
print(cluster_to_subtype)

lym.obs["lymphoid_subtype"] = (
    lym.obs["lymphoid_leiden"]
    .astype(str)
    .map(cluster_to_subtype)
)


# ============================================================
# 5. Transfer refined annotation to whole dataset
# ============================================================

adata.obs["cell_type_refined"] = (
    adata.obs["cell_type"]
    .astype(str)
)

adata.obs.loc[
    lym.obs_names,
    "cell_type_refined"
] = lym.obs["lymphoid_subtype"]

adata.obs["cell_type_refined"] = (
    adata.obs["cell_type_refined"]
    .astype("category")
)

print("\nRefined cell types:")
print(
    adata.obs["cell_type_refined"]
    .value_counts()
)


# ============================================================
# 6. Refined UMAP
# ============================================================

print("\n=== 6. Refined UMAP ===")

sc.pl.umap(
    adata,
    color="cell_type_refined",
    legend_loc="right margin",
    title="PBMC3k refined cell-type annotation",
    save="_day7_01_refined_celltypes.png",
)


# ============================================================
# 7. Differential marker analysis by refined cell type
# ============================================================

print("\n=== 7. Refined cell-type marker analysis ===")

sc.tl.rank_genes_groups(
    adata,
    groupby="cell_type_refined",
    method="wilcoxon",
    use_raw=True,
    key_added="rank_genes_refined",
)

marker_df = sc.get.rank_genes_groups_df(
    adata,
    group=None,
    key="rank_genes_refined",
)

marker_df.to_csv(
    OUT_DIR / "day7_refined_celltype_markers.csv",
    index=False,
)

print(
    marker_df[
        [
            "group",
            "names",
            "scores",
            "logfoldchanges",
            "pvals_adj",
        ]
    ].head(30)
)


# ============================================================
# 8. Canonical-marker DotPlot
# ============================================================

print("\n=== 8. Canonical marker DotPlot ===")

marker_dict = {

    "CD4 T": [
        "IL7R",
        "CCR7",
        "LTB",
    ],

    "CD8 T": [
        "CD3D",
        "CD8A",
        "CCL5",
        "GZMK",
    ],

    "NK": [
        "NKG7",
        "GNLY",
        "PRF1",
    ],

    "B": [
        "MS4A1",
        "CD79A",
        "CD74",
    ],

    "Monocyte": [
        "LYZ",
        "LST1",
        "S100A8",
        "S100A9",
    ],

    "Dendritic": [
        "FCER1A",
        "CST3",
        "HLA-DRA",
    ],

    "Platelet": [
        "PPBP",
        "PF4",
    ],

    "Cycling": [
        "TYMS",
        "PCNA",
    ],
}

available_markers = {}

for group, genes in marker_dict.items():

    valid = [
        g for g in genes
        if g in adata.raw.var_names
    ]

    if valid:
        available_markers[group] = valid

sc.pl.dotplot(
    adata,
    available_markers,
    groupby="cell_type_refined",
    use_raw=True,
    standard_scale="var",
    save="_day7_02_marker_dotplot.png",
)


# ============================================================
# 9. Marker heatmap
# ============================================================

print("\n=== 9. Marker heatmap ===")

top_marker_genes = []

for group in marker_df["group"].unique():

    sub = (
        marker_df[
            marker_df["group"] == group
        ]
        .sort_values(
            "scores",
            ascending=False,
        )
    )

    genes = (
        sub["names"]
        .dropna()
        .head(3)
        .tolist()
    )

    top_marker_genes.extend(genes)

# Remove duplicates while preserving order
top_marker_genes = list(
    dict.fromkeys(top_marker_genes)
)

sc.pl.heatmap(
    adata,
    top_marker_genes,
    groupby="cell_type_refined",
    use_raw=True,
    standard_scale="var",
    show_gene_labels=True,
    save="_day7_03_top_markers.png",
)


# ============================================================
# 10. Cell-type composition
# ============================================================

print("\n=== 10. Cell-type composition ===")

composition = (
    adata.obs["cell_type_refined"]
    .value_counts()
    .rename_axis("cell_type")
    .reset_index(name="n_cells")
)

composition["fraction"] = (
    composition["n_cells"]
    / composition["n_cells"].sum()
)

composition.to_csv(
    OUT_DIR / "day7_celltype_composition.csv",
    index=False,
)

print(composition)

fig, ax = plt.subplots(
    figsize=(8, 5)
)

ax.bar(
    composition["cell_type"],
    composition["fraction"],
)

ax.set_ylabel(
    "Fraction of cells"
)

ax.set_title(
    "PBMC3k cell-type composition"
)

ax.tick_params(
    axis="x",
    rotation=35,
)

fig.tight_layout()

fig.savefig(
    FIG_DIR / "day7_04_celltype_composition.png",
    dpi=200,
)

plt.close(fig)


# ============================================================
# 11. NK-like vs CD8 T-like differential expression
# ============================================================

print("\n=== 11. NK-like vs CD8 T-like DE ===")

contrast_groups = {
    "NK-like",
    "CD8 T-like",
}

present_groups = set(
    adata.obs["cell_type_refined"]
    .astype(str)
    .unique()
)

if contrast_groups.issubset(
    present_groups
):

    contrast = adata[
        adata.obs["cell_type_refined"]
        .isin(
            [
                "NK-like",
                "CD8 T-like",
            ]
        )
    ].copy()

    sc.tl.rank_genes_groups(
        contrast,
        groupby="cell_type_refined",
        groups=["NK-like"],
        reference="CD8 T-like",
        method="wilcoxon",
        use_raw=True,
        key_added="nk_vs_cd8",
    )

    nk_vs_cd8 = (
        sc.get.rank_genes_groups_df(
            contrast,
            group="NK-like",
            key="nk_vs_cd8",
        )
    )

    nk_vs_cd8.to_csv(
        OUT_DIR / "day7_NK_vs_CD8_DE.csv",
        index=False,
    )

    print(
        nk_vs_cd8[
            [
                "names",
                "scores",
                "logfoldchanges",
                "pvals_adj",
            ]
        ].head(20)
    )

else:

    print(
        "NK-like and CD8 T-like groups were not "
        "both detected; skipping this contrast."
    )


# ============================================================
# 12. Save
# ============================================================

print("\n=== 12. Save interpreted AnnData ===")

adata.write(
    OUT_DIR / "pbmc3k_day7_interpreted.h5ad"
)

adata.obs.to_csv(
    OUT_DIR / "day7_cell_metadata.csv"
)

cluster_scores.to_csv(
    OUT_DIR / "day7_lymphoid_signature_scores.csv"
)

print(
    "\nSaved:",
    OUT_DIR / "pbmc3k_day7_interpreted.h5ad"
)

print("\nDay 7 core analysis complete.")
