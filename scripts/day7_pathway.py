from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp


OUT = Path("outputs")
FIG = Path("figures")

OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)


# ============================================================
# 1. Load NK vs CD8 differential-expression results
# ============================================================

print("\n=== 1. Load NK vs CD8 DE ===")

de = pd.read_csv(
    OUT / "day7_NK_vs_CD8_DE.csv"
)

print(de.head())


# ============================================================
# 2. Define significant NK-high / CD8-high genes
# ============================================================

print("\n=== 2. Select genes ===")

nk_high = (
    de[
        (de["pvals_adj"] < 0.05)
        & (de["logfoldchanges"] > 0.5)
    ]
    .sort_values(
        "logfoldchanges",
        ascending=False,
    )
    ["names"]
    .dropna()
    .astype(str)
    .drop_duplicates()
    .head(200)
    .tolist()
)

cd8_high = (
    de[
        (de["pvals_adj"] < 0.05)
        & (de["logfoldchanges"] < -0.5)
    ]
    .sort_values(
        "logfoldchanges",
        ascending=True,
    )
    ["names"]
    .dropna()
    .astype(str)
    .drop_duplicates()
    .head(200)
    .tolist()
)

print("NK-high genes:", len(nk_high))
print("CD8-high genes:", len(cd8_high))

print("\nTop NK-high:")
print(nk_high[:20])

print("\nTop CD8-high:")
print(cd8_high[:20])


# ============================================================
# 3. Enrichr pathway databases
# ============================================================

gene_sets = [
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "KEGG_2021_Human",
]


def run_enrichment(label, genes):

    if len(genes) < 5:

        print(
            f"{label}: too few genes, skipping."
        )

        return None

    print(
        f"\n=== Enrichment: {label} ==="
    )

    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets,
        organism="human",
        outdir=None,
        cutoff=0.05,
    )

    res = enr.results.copy()

    res.to_csv(
        OUT / f"day7_{label}_pathway_enrichment.csv",
        index=False,
    )

    significant = (
        res[
            res["Adjusted P-value"] < 0.05
        ]
        .sort_values(
            "Adjusted P-value"
        )
    )

    print(
        significant[
            [
                "Gene_set",
                "Term",
                "Adjusted P-value",
                "Odds Ratio",
                "Combined Score",
            ]
        ]
        .head(20)
        .to_string(index=False)
    )

    return significant


nk_res = run_enrichment(
    "NK_high",
    nk_high,
)

cd8_res = run_enrichment(
    "CD8_high",
    cd8_high,
)


# ============================================================
# 4. Simple pathway plots
# ============================================================

def plot_top_pathways(
    result,
    label,
):

    if result is None or result.empty:
        return

    top = (
        result
        .head(15)
        .copy()
    )

    top["minus_log10_FDR"] = (
        -np.log10(
            top["Adjusted P-value"]
            .clip(lower=1e-300)
        )
    )

    top = top.iloc[::-1]

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    ax.barh(
        top["Term"],
        top["minus_log10_FDR"],
    )

    ax.set_xlabel(
        "-log10 adjusted P-value"
    )

    ax.set_title(
        f"{label} pathway enrichment"
    )

    fig.tight_layout()

    fig.savefig(
        FIG / f"day7_{label}_pathway_enrichment.png",
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)


plot_top_pathways(
    nk_res,
    "NK_high",
)

plot_top_pathways(
    cd8_res,
    "CD8_high",
)


print(
    "\nDay 7 pathway enrichment complete."
)
