# Build a Variance-Weighted Nearest-Neighbour Graph

Scales each principal component by a variance-derived weight (the same
weighting as
[`RunWeightedUMAP`](https://lachland.github.io/weightedUMAP/reference/RunWeightedUMAP.md)),
stores the result as a new `DimReduc`, and builds a KNN/SNN graph via
Seurat's
[`FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html).
The resulting graphs can be passed directly to
[`FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)
and, via the `graph` argument of
[`RunWeightedUMAP`](https://lachland.github.io/weightedUMAP/reference/RunWeightedUMAP.md),
for UMAP — ensuring that clustering and visualisation share exactly the
same nearest-neighbour structure.

## Usage

``` r
RunWeightedNeighbors(
  object,
  reduction = "pca",
  dims = NULL,
  weight.by = c("prop.var", "stdev", "mp", "none"),
  k.param = 20L,
  reduction.name = "wt.pca",
  reduction.key = "wtPCA_",
  prefix = "wt",
  weight.factor = 1,
  log.scale = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object with a PCA (or other linear) reduction already
  computed.

- reduction:

  Name of the dimensional reduction to use as input. Default: `"pca"`.

- dims:

  Integer vector of dimensions to use. If `NULL` (default), all
  available dimensions are used and a warning is raised when there are
  more than 50.

- weight.by:

  Weighting scheme applied to PC scores before UMAP. One of:

  `"prop.var"`

  :   Proportion of variance explained (`sdev^2 / sum(sdev^2)`).
      Default.

  `"stdev"`

  :   Standard deviation, normalised (`sdev / sum(sdev)`).

  `"none"`

  :   No weighting; equivalent to standard UMAP on PCA scores.

- k.param:

  Number of nearest neighbours for the KNN graph. Default: `20L`.

- reduction.name:

  Name under which to store the intermediate weighted-PC embedding in
  `object`. Default: `"wt.pca"`.

- reduction.key:

  Column-name prefix for the weighted-PC dimensions. Default:
  `"wtPCA_"`.

- prefix:

  Prefix for the two graphs stored in the Seurat object: `<prefix>_nn`
  (KNN, used for UMAP) and `<prefix>_snn` (SNN, used for clustering).
  Pass this prefix to `FindClusters(graph.name = "<prefix>_snn")` and
  `RunWeightedUMAP(graph = "<prefix>_nn")`. Default: `"wt"`.

- weight.factor:

  A number in `[0, 1]` controlling the influence of the weighting
  scheme. At `1` (default) full variance-derived weights are applied. At
  `0` all PCs contribute equally (standard nearest-neighbour graph).
  Must match the `weight.factor` used in the subsequent
  [`RunWeightedUMAP`](https://lachland.github.io/weightedUMAP/reference/RunWeightedUMAP.md)
  call to keep clustering and visualisation in the same space.

- log.scale:

  Logical. If `TRUE`, [`log1p()`](https://rdrr.io/r/base/Log.html) is
  applied to the weights before the `weight.factor` interpolation.
  Should match the value used in the subsequent
  [`RunWeightedUMAP`](https://lachland.github.io/weightedUMAP/reference/RunWeightedUMAP.md)
  call. Default: `FALSE`.

- verbose:

  Print progress messages and per-PC weights? Default: `TRUE`.

- ...:

  Additional arguments forwarded to
  [`FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html).

## Value

The input Seurat `object` with:

- A new `DimReduc` stored under `reduction.name` containing the
  variance-weighted PC embeddings.

- Two neighbour graphs (`<prefix>_nn` and `<prefix>_snn`) suitable for
  downstream
  [`FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)
  and `RunWeightedUMAP(graph = ...)`.

## See also

[`RunWeightedUMAP`](https://lachland.github.io/weightedUMAP/reference/RunWeightedUMAP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(wUMAP)

# 1. Compute the weighted KNN/SNN graphs (same space for clustering + UMAP)
pbmc <- RunWeightedNeighbors(pbmc, dims = 1:30, weight.by = "prop.var",
                             prefix = "wt")

# 2. Cluster on the weighted SNN graph
pbmc <- FindClusters(pbmc, graph.name = "wt_snn")

# 3. UMAP from the same weighted KNN graph — topology is identical to step 2
pbmc <- RunWeightedUMAP(pbmc, graph = "wt_nn", n.neighbors = 20,
                        reduction.name = "wt.umap")

DimPlot(pbmc, reduction = "wt.umap", label = TRUE)

# Use log.scale to compress the weight dynamic range so intermediate PCs
# contribute more alongside the dominant early PCs
pbmc <- RunWeightedNeighbors(pbmc, dims = 1:30, weight.by = "prop.var",
                             log.scale = TRUE, prefix = "wt.log")
pbmc <- FindClusters(pbmc, graph.name = "wt.log_snn")
pbmc <- RunWeightedUMAP(pbmc, graph = "wt.log_nn",
                        reduction.name = "wt.umap.log")
} # }
```
