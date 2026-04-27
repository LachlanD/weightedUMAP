# wUMAP: Variance-Weighted UMAP for Seurat Objects

`wUMAP` provides `RunWeightedUMAP()`, a drop-in replacement for Seurat's `RunUMAP()` that scales each principal component by a weight derived from its variance explained before running manifold learning. This gives biologically informative PCs more influence over the final embedding.

---

## Motivation

Standard UMAP treats every PC dimension equally when computing cell–cell distances. In practice, PC1 might explain 15% of variance while PC30 explains 0.3% — yet both contribute equally to the distance matrix. Weighting PCs by their variance explained amplifies axes of real biological variation and suppresses noise-dominated dimensions, often producing tighter, better-separated clusters.

---

## Installation

```r
# Install from GitHub (requires remotes)
remotes::install_github("your-username/weightedUMAP")
```

Dependencies (`uwot`, `SeuratObject`) will be installed automatically.

---

## Usage

### Basic usage

```r
library(Seurat)
library(wUMAP)

# Assumes pbmc has PCA already computed (e.g. via RunPCA())
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "pct.var")
DimPlot(pbmc, reduction = "wt.umap", label = TRUE)
```

### Weighting schemes

| `weight.by`   | Formula                              | Description                          |
|---------------|--------------------------------------|--------------------------------------|
| `"pct.var"`   | `sdev² / sum(sdev²) × 100`          | Percentage of variance explained (default) |
| `"prop.var"`  | `sdev² / sum(sdev²)`                | Proportion of variance explained     |
| `"eigenvalue"`| `sdev²`                             | Raw eigenvalue                       |
| `"stdev"`     | `sdev`                              | Standard deviation                   |
| `"none"`      | `1` (all equal)                      | No weighting — equivalent to standard UMAP |

```r
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "pct.var")    # default
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "stdev")
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "eigenvalue")
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "prop.var")
```

---

## Standard vs Weighted UMAP

The key difference is how PC scores are scaled before distances are computed:

- **Standard UMAP** (`weight.by = "none"`) — every PC contributes equally, so noisy high-index PCs can blur cluster boundaries.
- **Weighted UMAP** (`weight.by = "pct.var"`) — PCs are multiplied by their weight before UMAP runs, so dominant axes of biological variation drive the layout.

```r
library(patchwork)

# Run both on the same object
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "none",
                        reduction.name = "umap.std")
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "pct.var",
                        reduction.name = "wt.umap")

# Compare side by side
p1 <- DimPlot(pbmc, reduction = "umap.std", label = TRUE) +
  ggtitle("Standard UMAP")
p2 <- DimPlot(pbmc, reduction = "wt.umap",  label = TRUE) +
  ggtitle("Weighted UMAP (pct.var)")
p1 | p2
```

---

## Full argument reference

```r
RunWeightedUMAP(
  object,
  reduction      = "pca",       # source reduction
  dims           = NULL,        # e.g. 1:30; NULL uses all available dims
  weight.by      = "pct.var",   # weighting scheme (see table above)
  reduction.name = "wt.umap",   # name stored in the Seurat object
  reduction.key  = "wtUMAP_",   # column prefix: wtUMAP_1, wtUMAP_2, ...
  n.neighbors    = 30L,
  n.components   = 2L,
  metric         = "euclidean",
  min.dist       = 0.3,
  spread         = 1,
  seed.use       = 42L,
  verbose        = TRUE,
  ...                           # passed to uwot::umap()
)
```

The returned object gains a new `DimReduc` under `reduction.name`. Its `misc` slot stores the chosen `weight.by`, the `weights` vector, `source.reduction`, and `dims.used` for reproducibility.

---

## License

MIT
