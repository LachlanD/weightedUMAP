# Build a Local-PCA Nearest-Neighbour Graph

Computes per-cell local PCA distances (see
[`RunLocalPCAUMAP`](https://lachland.github.io/weightedUMAP/reference/RunLocalPCAUMAP.md))
and stores the resulting KNN and SNN graphs in the Seurat object. The
graphs can be passed to
[`FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)
(via the `graph.name` argument) and to
[`RunLocalPCAUMAP`](https://lachland.github.io/weightedUMAP/reference/RunLocalPCAUMAP.md)
(via the `nn_method` argument in `...`), ensuring that clustering and
visualisation share exactly the same local-PCA nearest-neighbour
structure.

## Usage

``` r
RunLocalPCANeighbors(
  object,
  reduction = "pca",
  dims = NULL,
  k.param = 20L,
  local.dims = NULL,
  local.weight.by = c("stdev", "prop.var", "none"),
  mp.filter = FALSE,
  prune.SNN = 1/15,
  prefix = "lp",
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object with PCA already run.

- reduction:

  Name of the PCA reduction to use. Default: `"pca"`.

- dims:

  Integer vector of PC dimensions to use. Default: all available.

- k.param:

  Number of nearest neighbours. Default: `20L`.

- local.dims:

  Number of local PCA directions to retain when computing distances.
  Must be \\\leq k.param\\. Default: `NULL` (uses `k.param` directions).

- local.weight.by:

  Weighting scheme for local PC directions. See
  [`RunLocalPCAUMAP`](https://lachland.github.io/weightedUMAP/reference/RunLocalPCAUMAP.md)
  for full descriptions. One of `"stdev"` (default), `"prop.var"`,
  `"none"`.

- mp.filter:

  Logical. If `TRUE`, local PC directions at or below the
  per-neighbourhood Marchenko-Pastur noise threshold are zeroed out
  after the base weights are computed. Can be combined with any
  `local.weight.by` scheme. Default: `FALSE`.

- prune.SNN:

  Jaccard similarity threshold below which edges are removed from the
  SNN graph. Default: `1/15`.

- prefix:

  Prefix for the two stored graphs: `<prefix>_nn` (KNN, binary) and
  `<prefix>_snn` (SNN, Jaccard-weighted). Pass `<prefix>_snn` to
  `FindClusters(graph.name = ...)`. Default: `"lp"`.

- verbose:

  Print progress messages. Default: `TRUE`.

## Value

The input Seurat `object` with two new Graph objects:

- `<prefix>_nn` — binary KNN graph (symmetric).

- `<prefix>_snn` — Jaccard-weighted SNN graph (sparse, pruned).

## Details

The SNN graph is a sparse matrix of Jaccard similarities, pruned at
`prune.SNN` (default \\1/15\\), matching Seurat's convention.

## See also

[`RunLocalPCAUMAP`](https://lachland.github.io/weightedUMAP/reference/RunLocalPCAUMAP.md),
[`RunWeightedNeighbors`](https://lachland.github.io/weightedUMAP/reference/RunWeightedNeighbors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(wUMAP)

# Assumes pbmc has PCA already run

# 1. Build local PCA KNN/SNN graphs — stdev weighting (default, recommended)
pbmc <- RunLocalPCANeighbors(pbmc, dims = 1:30, k.param = 20,
                             prefix = "lp")

# 2. Cluster on the local PCA SNN graph
pbmc <- FindClusters(pbmc, graph.name = "lp_snn", resolution = 0.5)

# 3. Visualise with matching local PCA UMAP (same k, same weighting)
pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 20,
                        reduction.name = "lp.umap")

DimPlot(pbmc, reduction = "lp.umap", label = TRUE)

# With MP filtering: discard local noise directions, then weight by stdev
pbmc <- RunLocalPCANeighbors(pbmc, dims = 1:30, k.param = 20,
                             mp.filter = TRUE, prefix = "lp.mp")
pbmc <- FindClusters(pbmc, graph.name = "lp.mp_snn", resolution = 0.5)
pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 20,
                        mp.filter = TRUE, reduction.name = "lp.umap.mp")
} # }
```
