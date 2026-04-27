# wUMAP: Variance-Weighted UMAP for Seurat Objects

`wUMAP` provides `RunWeightedUMAP()`, a drop-in replacement for Seurat's `RunUMAP()` that scales each principal component by a weight derived from its variance explained before running manifold learning. This gives biologically informative PCs more influence over the final embedding.

---

## Motivation

Standard UMAP treats every PC dimension equally when computing cell–cell distances. In practice, PC1 might explain 15% of variance while PC30 explains 0.3% — yet both contribute equally to the distance matrix. Weighting PCs by their variance explained amplifies axes of real biological variation and suppresses noise-dominated dimensions.

---

## Installation

```r
# Install from GitHub (requires remotes)
remotes::install_github("LachlanD/weightedUMAP")
```

Dependencies (`uwot`, `SeuratObject`) will be installed automatically.

---

## Usage

### Basic usage

```r
library(Seurat)
library(wUMAP)

# Assumes pbmc has PCA already computed (e.g. via RunPCA())
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "prop.var")
DimPlot(pbmc, reduction = "wt.umap", label = TRUE)
```

### `weight.factor`: continuous control

`weight.factor` (0–1) blends between standard UMAP (0) and fully weighted UMAP (1):

```r
# Standard UMAP — all PCs equal
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "prop.var",
                        weight.factor = 0)

# Half-strength weighting
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "prop.var",
                        weight.factor = 0.5)

# Full weighting (default)
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "prop.var",
                        weight.factor = 1)
```

### Weighting schemes

| `weight.by`   | Weight applied to PC _i_              | Description                          |
|---------------|---------------------------------------|--------------------------------------|
| `"prop.var"`  | `sdevᵢ² / Σ sdev²`                   | Proportion of variance explained (default) |
| `"stdev"`     | `sdevᵢ`                               | Standard deviation                   |
| `"none"`      | `1`                                   | No weighting — equivalent to standard UMAP |

```r
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "prop.var")    # default
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "stdev")
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "none")
```

---

## Standard vs Weighted UMAP

The key difference is how PC scores are scaled before distances are computed:

- **Standard UMAP** (`weight.by = "none"`) — every PC contributes equally, so noisy high-index PCs can blur cluster boundaries.
- **Weighted UMAP** (`weight.by = "prop.var"`) — PCs are multiplied by their weight before UMAP runs, so dominant axes of biological variation drive the layout.

```r
library(patchwork)

# Run both on the same object
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "none",
                        reduction.name = "umap.std")
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "prop.var",
                        reduction.name = "wt.umap")

# Compare side by side
p1 <- DimPlot(pbmc, reduction = "umap.std", label = TRUE) +
  ggtitle("Standard UMAP")
p2 <- DimPlot(pbmc, reduction = "wt.umap",  label = TRUE) +
  ggtitle("Weighted UMAP (prop.var)")
p1 | p2
```

---

## Full argument reference

```r
RunWeightedUMAP(
  object,
  reduction      = "pca",       # source reduction
  dims           = NULL,        # e.g. 1:30; NULL uses all available dims
  weight.by      = "prop.var",   # weighting scheme (see table above)
  weight.factor  = 1,           # 0 = standard UMAP, 1 = fully weighted
  graph          = NULL,        # name of KNN graph from RunWeightedNeighbors()
                                #   e.g. "wt_nn"; overrides dims/weight.by/reduction
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

## Consistent clustering and UMAP

By default, `FindNeighbors()` and `RunUMAP()` each build their own neighbour
graph — on unweighted PCA scores — so the topology used for clustering and the
one used for the UMAP embedding can differ.

`RunWeightedNeighbors()` computes the weighted PC embeddings **once**, stores
them as a `"<prefix>.pca"` reduction, and builds KNN/SNN graphs in that same
weighted space. When you call `RunWeightedUMAP(graph = "wt_nn", ...)`, it
looks up that stored weighted embedding and runs UMAP on it directly — so
clustering and visualisation both operate in the same weighted PC space with
the same `k`.

> **Note:** the `graph` argument tells `RunWeightedUMAP` *which weighted space*
> to use (by convention `"wt_nn"` → `"wt.pca"`). The actual UMAP is computed
> from the embedding, not from the binary adjacency values in the graph.

```r
library(Seurat)
library(wUMAP)
library(patchwork)

# Step 1 — build weighted KNN/SNN graphs (k = 20 by default)
#           stores 'wt.pca' embedding and 'wt_nn' / 'wt_snn' graphs
pbmc <- RunWeightedNeighbors(pbmc, dims = 1:30, weight.by = "prop.var",
                              k.param = 20, graph.name = "wt")

# Step 2 — cluster on the weighted SNN graph
pbmc <- FindClusters(pbmc, graph.name = "wt_snn")

# Step 3 — UMAP in the same weighted space
#           set n.neighbors to match the k.param used above
pbmc <- RunWeightedUMAP(pbmc, graph = "wt_nn", n.neighbors = 20,
                         reduction.name = "wt.umap")

DimPlot(pbmc, reduction = "wt.umap", label = TRUE)
```

Compare with the standard (inconsistent) workflow:

```r
# Standard workflow — neighbour graphs computed independently
pbmc <- FindNeighbors(pbmc, dims = 1:30)           # unweighted, k = 20
pbmc <- FindClusters(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:30)                 # separate NN computation
```

---

## License

MIT
