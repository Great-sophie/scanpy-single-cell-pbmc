from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp


# ============================================================
# Configuration
# ============================================================

OUT = Path("outputs")
FIG = Path("figures")

OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

DE_FILE = OUT / "day7_NK_vs_CD8_DE.csv"

GENE_SETS = [
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "KEGG_2021_Human",
]


# ============================================================
# 1. Load differential-expression results
# ============================================================

print("\n=== 1. Load NK vs CD8 differential expression ===")

if not DE_FILE.exists():
    raise FileNotFoundError(
        f"Cannot find: {DE_FILE}\n"
        "Run day7_refine_cytotoxic.py first."
    )

de = pd.read_csv(DE_FILE)

required_columns = {
    "names",
    "logfoldchanges",
    "pvals_adj",
}

missing = required_columns - set(de.columns)

if missing:
    raise ValueError(
        f"DE table missing columns: {sorted(missing)}"
    )

print("DE genes:", len(de))

print(
    de[
        [
            "names",
            "logfoldchanges",
            "pvals_adj",
        ]
    ]
    .head(10)
    .to_string(index=False)
)


# ============================================================
# 2. Housekeeping / technical gene filter
# ============================================================

print("\n=== 2. Define technical-gene filter ===")


def is_biological_gene(gene):
    """
    Remove genes that commonly dominate enrichment because of
    ribosomal, mitochondrial, histone, or generic technical programs.

    This is an exploratory filtering step for biological interpretation.
    The original unfiltered DE results remain unchanged.
    """

    gene = str(gene).strip().upper()

    excluded_prefixes = (
        "RPL",      # ribosomal protein large subunit
        "RPS",      # ribosomal protein small subunit
        "MT-",      # mitochondrial genes
        "HIST",     # histone genes
        "H1-",      # histone H1 family
    )

    if gene.startswith(excluded_prefixes):
        return False

    excluded_exact = {
        "MALAT1",
        "TPT1",
        "EEF1A1",
        "EEF1B2",
        "EEF2",
        "B2M",
        "GAPDH",
        "ACTB",
        "ACTG1",
        "PPIA",
        "FTL",
        "FTH1",
    }

    if gene in excluded_exact:
        return False

    return True


de["keep_for_filtered_pathway"] = (
    de["names"]
    .astype(str)
    .apply(is_biological_gene)
)

print(
    "Genes retained after filter:",
    int(de["keep_for_filtered_pathway"].sum()),
    "/",
    len(de),
)


# Save annotated DE table
de.to_csv(
    OUT / "day7_NK_vs_CD8_DE_with_pathway_filter.csv",
    index=False,
)


# ============================================================
# 3. Select significant NK-high genes
# ============================================================

print("\n=== 3. Select NK-high genes ===")

nk_table = (
    de[
        (de["pvals_adj"] < 0.05)
        & (de["logfoldchanges"] > 0.5)
        & (de["keep_for_filtered_pathway"])
    ]
    .sort_values(
        "logfoldchanges",
        ascending=False,
    )
    .copy()
)

nk_high = (
    nk_table["names"]
    .dropna()
    .astype(str)
    .drop_duplicates()
    .head(200)
    .tolist()
)

print(
    "NK-high genes retained:",
    len(nk_high),
)

print("\nTop NK-high genes:")

print(
    nk_table[
        [
            "names",
            "logfoldchanges",
            "pvals_adj",
        ]
    ]
    .head(30)
    .to_string(index=False)
)


# ============================================================
# 4. Select significant CD8-high genes
# ============================================================

print("\n=== 4. Select CD8-high genes ===")

cd8_table = (
    de[
        (de["pvals_adj"] < 0.05)
        & (de["logfoldchanges"] < -0.5)
        & (de["keep_for_filtered_pathway"])
    ]
    .sort_values(
        "logfoldchanges",
        ascending=True,
    )
    .copy()
)

cd8_high = (
    cd8_table["names"]
    .dropna()
    .astype(str)
    .drop_duplicates()
    .head(200)
    .tolist()
)

print(
    "CD8-high genes retained:",
    len(cd8_high),
)

print("\nTop CD8-high genes:")

print(
    cd8_table[
        [
            "names",
            "logfoldchanges",
            "pvals_adj",
        ]
    ]
    .head(30)
    .to_string(index=False)
)


# ============================================================
# 5. Save filtered gene lists
# ============================================================

print("\n=== 5. Save filtered gene lists ===")

pd.DataFrame(
    {
        "gene": nk_high
    }
).to_csv(
    OUT / "day7_NK_high_filtered_genes.csv",
    index=False,
)

pd.DataFrame(
    {
        "gene": cd8_high
    }
).to_csv(
    OUT / "day7_CD8_high_filtered_genes.csv",
    index=False,
)


# ============================================================
# 6. Enrichment function
# ============================================================

def run_enrichment(
    label,
    genes,
):
    """
    Run Enrichr over GO BP, Reactome and KEGG.
    """

    if len(genes) < 5:
        print(
            f"\n{label}: fewer than 5 genes; "
            "skipping enrichment."
        )
        return None

    print(
        f"\n=== 6. Enrichment: {label} ==="
    )

    print(
        "Genes submitted:",
        len(genes),
    )

    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=GENE_SETS,
        organism="human",
        outdir=None,
        cutoff=1.0,
    )

    result = enr.results.copy()

    if result.empty:
        print(
            f"No enrichment results returned for {label}."
        )
        return result

    result = result.sort_values(
        "Adjusted P-value",
        ascending=True,
    )

    result.to_csv(
        OUT
        / f"day7_{label}_filtered_pathway_enrichment.csv",
        index=False,
    )

    significant = result[
        result["Adjusted P-value"] < 0.05
    ].copy()

    print(
        f"Significant pathways (FDR < 0.05): "
        f"{len(significant)}"
    )

    if significant.empty:
        print(
            "No pathways reached FDR < 0.05."
        )
    else:

        cols = [
            "Gene_set",
            "Term",
            "Adjusted P-value",
            "Odds Ratio",
            "Combined Score",
        ]

        cols = [
            c
            for c in cols
            if c in significant.columns
        ]

        print(
            "\nTop 20 significant pathways:"
        )

        print(
            significant[cols]
            .head(20)
            .to_string(index=False)
        )

    return result


# ============================================================
# 7. Run enrichment
# ============================================================

nk_result = run_enrichment(
    "NK_high",
    nk_high,
)

cd8_result = run_enrichment(
    "CD8_high",
    cd8_high,
)


# ============================================================
# 8. Filter biologically informative pathways for plotting
# ============================================================

def clean_pathway_results(result):
    """
    Remove some highly generic terms only for visualization.
    The complete enrichment CSV is always preserved.
    """

    if result is None or result.empty:
        return result

    x = result.copy()

    x = x[
        x["Adjusted P-value"] < 0.05
    ].copy()

    if x.empty:
        return x

    generic_terms = [
        "immune system",
        "adaptive immune system",
        "innate immune system",
        "disease",
        "infection",
    ]

    # Keep broad terms in CSV, but deprioritize/remove the
    # most generic ones in the summary plot.
    keep = np.ones(
        len(x),
        dtype=bool,
    )

    terms_lower = (
        x["Term"]
        .astype(str)
        .str.lower()
    )

    for phrase in generic_terms:
        keep &= ~(
            terms_lower
            == phrase.lower()
        )

    cleaned = x.loc[keep].copy()

    # If filtering is too aggressive, fall back to all significant terms.
    if len(cleaned) < 5:
        cleaned = x.copy()

    return cleaned


# ============================================================
# 9. Plot pathway enrichment
# ============================================================

def plot_pathways(
    result,
    label,
):
    if result is None or result.empty:
        print(
            f"{label}: no result available for plotting."
        )
        return

    plot_df = clean_pathway_results(
        result
    )

    if plot_df is None or plot_df.empty:
        print(
            f"{label}: no significant pathways to plot."
        )
        return

    top = (
        plot_df
        .sort_values(
            "Adjusted P-value",
            ascending=True,
        )
        .head(15)
        .copy()
    )

    top["minus_log10_FDR"] = (
        -np.log10(
            top["Adjusted P-value"]
            .clip(lower=1e-300)
        )
    )

    # Reverse order so strongest pathway appears at top
    top = top.iloc[::-1]

    fig, ax = plt.subplots(
        figsize=(10, 7)
    )

    ax.barh(
        top["Term"],
        top["minus_log10_FDR"],
    )

    ax.set_xlabel(
        "-log10 adjusted P-value"
    )

    ax.set_ylabel(
        ""
    )

    ax.set_title(
        f"{label}: filtered pathway enrichment"
    )

    fig.tight_layout()

    output_path = (
        FIG
        / f"day7_{label}_filtered_pathway_enrichment.png"
    )

    fig.savefig(
        output_path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(
        "Saved figure:",
        output_path,
    )


plot_pathways(
    nk_result,
    "NK_high",
)

plot_pathways(
    cd8_result,
    "CD8_high",
)


# ============================================================
# 10. Make compact summary tables
# ============================================================

def save_top_summary(
    result,
    label,
):
    if result is None or result.empty:
        return

    sig = (
        result[
            result["Adjusted P-value"] < 0.05
        ]
        .sort_values(
            "Adjusted P-value"
        )
        .copy()
    )

    if sig.empty:
        return

    cols = [
        "Gene_set",
        "Term",
        "Adjusted P-value",
        "Odds Ratio",
        "Combined Score",
        "Overlap",
        "Genes",
    ]

    cols = [
        c
        for c in cols
        if c in sig.columns
    ]

    sig[
        cols
    ].head(50).to_csv(
        OUT
        / f"day7_{label}_filtered_pathway_TOP50.csv",
        index=False,
    )


save_top_summary(
    nk_result,
    "NK_high",
)

save_top_summary(
    cd8_result,
    "CD8_high",
)


# ============================================================
# 11. Final summary
# ============================================================

print(
    "\n========================================"
)

print(
    "Day 7 filtered pathway enrichment complete."
)

print(
    "========================================"
)

print("\nGenerated outputs:")

print(
    OUT
    / "day7_NK_high_filtered_pathway_enrichment.csv"
)

print(
    OUT
    / "day7_CD8_high_filtered_pathway_enrichment.csv"
)

print(
    OUT
    / "day7_NK_high_filtered_pathway_TOP50.csv"
)

print(
    OUT
    / "day7_CD8_high_filtered_pathway_TOP50.csv"
)

print("\nGenerated figures:")

print(
    FIG
    / "day7_NK_high_filtered_pathway_enrichment.png"
)

print(
    FIG
    / "day7_CD8_high_filtered_pathway_enrichment.png"
)
