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
  distances. Must be \\\leq\\ `k.param`. Default: `NULL` (uses all
  `k.param` directions).

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
`reduction.name`.

## Algorithm

1.  Find the global `k.param` nearest neighbours of every cell in the
    supplied PCA embedding (using RANN).

2.  For each cell \\i\\, collect its \\k\\-neighbourhood, centre by the
    neighbourhood mean, and compute a compact SVD.

3.  Re-express the displacement vectors from \\i\\ to each neighbour in
    the local `local.dims` principal directions.

4.  Report the Euclidean norm in local-PC space as the refined distance.

5.  Pass the \\n \times k\\ index/distance matrices directly to
    [`umap`](https://jlmelville.github.io/uwot/reference/umap.html) as a
    precomputed k-NN graph.

## See also

[`RunWeightedUMAP`](https://lachland.github.io/weightedUMAP/reference/RunWeightedUMAP.md),
[`RunWeightedNeighbors`](https://lachland.github.io/weightedUMAP/reference/RunWeightedNeighbors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(wUMAP)

# Assumes pbmc has PCA already run
pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 30,
                        reduction.name = "lp.umap")

# Compare to standard UMAP
p1 <- DimPlot(pbmc, reduction = "umap",    label = TRUE) + ggtitle("Standard")
p2 <- DimPlot(pbmc, reduction = "lp.umap", label = TRUE) + ggtitle("Local PCA")
p1 | p2
} # }
```
