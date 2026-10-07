from pathlib import Path

import scanpy as sc
import pandas as pd


OUT = Path("outputs")
FIG = Path("figures")

sc.settings.figdir = str(FIG)
sc.settings.verbosity = 2


# ============================================================
# 1. Load Day 6 annotated data
# ============================================================

print("\n=== 1. Load Day 6 object ===")

adata = sc.read_h5ad(
    OUT / "pbmc3k_annotated.h5ad"
)

print(adata)


# ============================================================
# 2. Extract mixed cytotoxic compartment
# ============================================================

print("\n=== 2. Extract NK / cytotoxic compartment ===")

mask = (
    adata.obs["cell_type"]
    .astype(str)
    == "NK / cytotoxic lymphocytes"
)

cyt_original = adata[mask].copy()

print("Cytotoxic compartment:", cyt_original.n_obs)


# Use full normalized/log1p expression stored in raw
cyt = cyt_original.raw.to_adata()

# preserve metadata
cyt.obs = cyt_original.obs.copy()

# preserve full expression for marker testing
cyt.raw = cyt.copy()


# ============================================================
# 3. Score T / CD8 / NK programs
# ============================================================

print("\n=== 3. Signature scoring ===")

signatures = {

    "T_identity": [
        "CD3D",
        "CD3E",
        "TRAC",
        "CD247",
        "LCK",
    ],

    "CD8_cytotoxic": [
        "CD8A",
        "CD8B",
        "CCL5",
        "GZMK",
        "CTSW",
    ],

    "NK_identity": [
        "NKG7",
        "GNLY",
        "KLRD1",
        "TYROBP",
        "FCER1G",
    ],
}


for name, genes in signatures.items():

    valid = [
        g for g in genes
        if g in cyt.var_names
    ]

    print(name, ":", valid)

    sc.tl.score_genes(
        cyt,
        gene_list=valid,
        score_name=f"{name}_score",
        random_state=42,
    )


# ============================================================
# 4. HVG
# ============================================================

print("\n=== 4. HVG ===")

sc.pp.highly_variable_genes(
    cyt,
    n_top_genes=1000,
    flavor="seurat",
)

print(
    "HVG:",
    int(cyt.var["highly_variable"].sum())
)

cyt = cyt[
    :,
    cyt.var["highly_variable"],
].copy()


# ============================================================
# 5. Scale + PCA
# ============================================================

print("\n=== 5. PCA ===")

sc.pp.scale(
    cyt,
    max_value=10,
)

sc.tl.pca(
    cyt,
    n_comps=30,
    svd_solver="arpack",
)


# ============================================================
# 6. Neighbor graph + UMAP
# ============================================================

print("\n=== 6. Neighbors / UMAP ===")

sc.pp.neighbors(
    cyt,
    n_neighbors=10,
    n_pcs=20,
)

sc.tl.umap(
    cyt,
    random_state=42,
)


# ============================================================
# 7. Higher-resolution Leiden
# ============================================================

print("\n=== 7. Cytotoxic subclustering ===")

sc.tl.leiden(
    cyt,
    resolution=0.8,
    key_added="cytotoxic_leiden",
    random_state=42,
)

print(
    cyt.obs["cytotoxic_leiden"]
    .value_counts()
    .sort_index()
)


# ============================================================
# 8. Cluster markers
# ============================================================

print("\n=== 8. Cytotoxic cluster markers ===")

sc.tl.rank_genes_groups(
    cyt,
    groupby="cytotoxic_leiden",
    method="wilcoxon",
    use_raw=True,
    key_added="cytotoxic_markers",
)

markers = sc.get.rank_genes_groups_df(
    cyt,
    group=None,
    key="cytotoxic_markers",
)

markers.to_csv(
    OUT / "day7_cytotoxic_subcluster_markers.csv",
    index=False,
)


for group in sorted(
    markers["group"].astype(str).unique(),
    key=int,
):

    sub = (
        markers[
            markers["group"].astype(str) == group
        ]
        .sort_values(
            "scores",
            ascending=False,
        )
    )

    print(
        f"\n--- Cytotoxic cluster {group} ---"
    )

    print(
        sub[
            [
                "names",
                "scores",
                "logfoldchanges",
                "pvals_adj",
            ]
        ]
        .head(10)
        .to_string(index=False)
    )


# ============================================================
# 9. Mean signature score per cluster
# ============================================================

score_cols = [
    "T_identity_score",
    "CD8_cytotoxic_score",
    "NK_identity_score",
]

cluster_scores = (
    cyt.obs
    .groupby(
        "cytotoxic_leiden",
        observed=True,
    )[score_cols]
    .mean()
)

print("\n=== Mean scores per cytotoxic cluster ===")

print(cluster_scores)

cluster_scores.to_csv(
    OUT / "day7_cytotoxic_cluster_scores.csv"
)


# ============================================================
# 10. Conservative CD8 vs NK annotation
# ============================================================

cluster_to_type = {}

for cluster, row in cluster_scores.iterrows():

    # Strong T lineage identity supports CD8 T
    # Strong NK identity supports NK
    if (
        row["T_identity_score"]
        > row["NK_identity_score"]
    ):

        label = "CD8 T-like"

    else:

        label = "NK-like"

    cluster_to_type[str(cluster)] = label


print("\n=== Cytotoxic cluster → subtype ===")

print(cluster_to_type)


cyt.obs["cytotoxic_subtype"] = (
    cyt.obs["cytotoxic_leiden"]
    .astype(str)
    .map(cluster_to_type)
)


print("\nSubtype counts:")

print(
    cyt.obs["cytotoxic_subtype"]
    .value_counts()
)


# ============================================================
# 11. Cytotoxic UMAP
# ============================================================

sc.pl.umap(
    cyt,
    color=[
        "cytotoxic_leiden",
        "cytotoxic_subtype",
        "T_identity_score",
        "CD8_cytotoxic_score",
        "NK_identity_score",
    ],
    ncols=2,
    save="_day7_05_cytotoxic_refinement.png",
)


# ============================================================
# 12. Transfer labels to complete dataset
# ============================================================

print("\n=== 12. Transfer refined labels ===")

adata.obs["cell_type_refined_v2"] = (
    adata.obs["cell_type"]
    .astype(str)
)

# Broad T cluster is kept conservative
adata.obs.loc[
    adata.obs["cell_type"]
    .astype(str)
    == "T cells",
    "cell_type_refined_v2",
] = "Conventional T cells"


# Assign cytotoxic cells
adata.obs.loc[
    cyt.obs_names,
    "cell_type_refined_v2",
] = cyt.obs[
    "cytotoxic_subtype"
].astype(str)


adata.obs["cell_type_refined_v2"] = (
    adata.obs["cell_type_refined_v2"]
    .astype("category")
)


print(
    adata.obs["cell_type_refined_v2"]
    .value_counts()
)


# ============================================================
# 13. Global UMAP
# ============================================================

sc.pl.umap(
    adata,
    color="cell_type_refined_v2",
    legend_loc="right margin",
    title="PBMC3k refined annotation",
    save="_day7_06_global_refined_v2.png",
)


# ============================================================
# 14. NK vs CD8 differential expression
# ============================================================

print("\n=== 14. NK vs CD8 differential expression ===")

present = set(
    adata.obs[
        "cell_type_refined_v2"
    ]
    .astype(str)
    .unique()
)

if {
    "NK-like",
    "CD8 T-like",
}.issubset(present):

    contrast = adata[
        adata.obs[
            "cell_type_refined_v2"
        ].isin(
            [
                "NK-like",
                "CD8 T-like",
            ]
        )
    ].copy()

    print(
        contrast.obs[
            "cell_type_refined_v2"
        ].value_counts()
    )

    sc.tl.rank_genes_groups(
        contrast,
        groupby="cell_type_refined_v2",
        groups=["NK-like"],
        reference="CD8 T-like",
        method="wilcoxon",
        use_raw=True,
        key_added="NK_vs_CD8",
    )

    de = sc.get.rank_genes_groups_df(
        contrast,
        group="NK-like",
        key="NK_vs_CD8",
    )

    de.to_csv(
        OUT / "day7_NK_vs_CD8_DE.csv",
        index=False,
    )

    print("\n=== Top NK-high genes ===")

    print(
        de[
            [
                "names",
                "scores",
                "logfoldchanges",
                "pvals_adj",
            ]
        ]
        .head(20)
        .to_string(index=False)
    )

else:

    print(
        "Still unable to separate NK-like "
        "and CD8 T-like populations."
    )


# ============================================================
# 15. Save
# ============================================================

adata.write(
    OUT / "pbmc3k_day7_refined_v2.h5ad"
)

cyt.write(
    OUT / "pbmc3k_cytotoxic_subcluster.h5ad"
)

print("\nSaved refined objects.")
print("\nDay 7 cytotoxic refinement complete.")
