# UMAP with Local PCA Distances

Computes UMAP using a **local PCA metric**: for each cell, the distance
to each of its \\k\\ nearest neighbours is re-measured in a
locally-fitted PCA basis built from that neighbourhood. This captures
anisotropic local structure (e.g. the tangent direction of a trajectory)
and de-emphasises transverse noise, producing a richer distance measure
than a single global Euclidean metric.

## Usage

``` r
RunLocalPCAUMAP(
  object,
  reduction = "pca",
  dims = NULL,
  k.param = 30L,
  local.dims = NULL,
  local.weight.by = c("stdev", "prop.var", "none"),
  mp.filter = FALSE,
  reduction.name = "lp.umap",
  reduction.key = "lpUMAP_",
  n.components = 2L,
  min.dist = 0.3,
  spread = 1,
  seed.use = 42L,
  verbose = TRUE,
  ...
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

  Number of nearest neighbours for the local PCA. Default: `30L`.

- local.dims:

  Number of local PCA directions to retain when computing refined
  distances. Must be \\\leq\\ `k.param`. Default: `NULL` (uses `k.param`
  directions).

- local.weight.by:

  Weighting scheme applied to local PC directions. Analogous to
  `weight.by` in
  [`RunWeightedUMAP`](https://lachland.github.io/weightedUMAP/reference/RunWeightedUMAP.md),
  but applied to the *local* singular values of each neighbourhood SVD.

  `"stdev"`

  :   Weight by local standard deviation (\\s_l / \sum s\\). Default.
      Gently emphasises the dominant local direction without collapsing
      minor variation.

  `"prop.var"`

  :   Weight each local PC by its proportion of local variance (\\s_l^2
      / \sum s^2\\); more aggressively emphasises the leading local
      direction.

  `"none"`

  :   No weighting; equal contribution from all local PCs.

- mp.filter:

  Logical. If `TRUE`, local PC directions whose local variance is at or
  below the per-neighbourhood Marchenko-Pastur noise threshold
  (\\\lambda\_{\max} = (1 + \sqrt{d/k})^2\\) are zeroed out after the
  `local.weight.by` weights are computed, and the remaining weights are
  renormalised. Can be combined with any `local.weight.by` scheme.
  Default: `FALSE`.

- reduction.name:

  Name under which to store the UMAP result. Default: `"lp.umap"`.

- reduction.key:

  Column-name prefix for the UMAP dimensions. Default: `"lpUMAP_"`.

- n.components:

  Number of UMAP dimensions. Default: `2L`.

- min.dist:

  UMAP `min_dist` parameter. Default: `0.3`.

- spread:

  UMAP `spread` parameter. Default: `1`.

- seed.use:

  Random seed. Default: `42L`.

- verbose:

  Print progress messages. Default: `TRUE`.

- ...:

  Additional arguments forwarded to
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.html).

## Value

The input Seurat `object` with a new DimReduc stored under
`reduction.name`. The `misc` slot records `source.reduction`,
`dims.used`, `k.param`, `local.dims`, and `local.weight.by`.

## Details

Optionally, each local PC direction can be weighted by its local
variance contribution (`local.weight.by`), analogous to how `weight.by`
weights global PCs in
[`RunWeightedUMAP`](https://lachland.github.io/weightedUMAP/reference/RunWeightedUMAP.md).

## Algorithm

1.  Find the global `k.param` nearest neighbours of every cell in the
    supplied PCA embedding (using RANN).

2.  For each cell \\i\\, collect its \\k\\-neighbourhood, centre by the
    neighbourhood mean, and compute a compact SVD.

3.  Optionally weight each of the `local.dims` local PC directions by
    its local variance contribution (`local.weight.by`).

4.  Report the (weighted) Euclidean norm in local-PC space as the
    refined distance.

5.  Pass the \\n \times k\\ index/distance matrices directly to
    [`umap`](https://jlmelville.github.io/uwot/reference/umap.html) as a
    precomputed k-NN graph.

## See also

[`RunLocalPCANeighbors`](https://lachland.github.io/weightedUMAP/reference/RunLocalPCANeighbors.md),
[`RunWeightedUMAP`](https://lachland.github.io/weightedUMAP/reference/RunWeightedUMAP.md),
[`RunWeightedNeighbors`](https://lachland.github.io/weightedUMAP/reference/RunWeightedNeighbors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(wUMAP)

# Assumes pbmc has PCA already run

# Local PCA UMAP with stdev weighting of local PC directions (default)
pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 30,
                        reduction.name = "lp.umap")

# Unweighted local PCA for comparison
pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 30,
                        local.weight.by = "none",
                        reduction.name = "lp.umap.unweighted")

# With MP filtering: zero local noise directions, then apply stdev weights
pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 30,
                        mp.filter = TRUE,
                        reduction.name = "lp.umap.mp")

# Compare standard, unweighted local PCA, and weighted local PCA
p1 <- DimPlot(pbmc, reduction = "umap",               label = TRUE) +
  ggtitle("Standard")
p2 <- DimPlot(pbmc, reduction = "lp.umap.unweighted", label = TRUE) +
  ggtitle("Local PCA (unweighted)")
p3 <- DimPlot(pbmc, reduction = "lp.umap",            label = TRUE) +
  ggtitle("Local PCA (stdev)")
p1 | p2 | p3
} # }
```
