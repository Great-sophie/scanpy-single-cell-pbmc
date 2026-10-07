# Scanpy Single-Cell RNA-seq Analysis — PBMC3k

End-to-end single-cell RNA-seq analysis using **Scanpy** and the PBMC3k dataset, covering quality control, normalization, dimensionality reduction, graph-based clustering, marker-based annotation, hierarchical subclustering, differential expression, and pathway enrichment.

## Overview

The analysis workflow is:

```text
Raw counts
   ↓
QC / filtering
   ↓
Normalization + log1p
   ↓
Highly variable genes
   ↓
PCA
   ↓
KNN graph
   ├── UMAP visualization
   └── Leiden clustering
          ↓
     Marker discovery
          ↓
   Cell-type annotation
          ↓
Hierarchical subclustering
          ↓
CD8 T-like vs NK-like refinement
          ↓
Differential expression
          ↓
Pathway enrichment
          ↓
Biological interpretation

