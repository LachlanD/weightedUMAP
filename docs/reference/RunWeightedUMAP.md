# Run Variance-Weighted UMAP on a Seurat Object

Scales each principal component embedding by a weight derived from the
variance explained by that component, then runs UMAP via uwot. The
result is stored as a new `DimReduc` in the returned Seurat object.

## Usage

``` r
RunWeightedUMAP(
  object,
  reduction = "pca",
  dims = NULL,
  weight.by = c("prop.var", "stdev", "none"),
  graph = NULL,
  reduction.name = "wt.umap",
  reduction.key = "wtUMAP_",
  n.neighbors = 30L,
  n.components = 2L,
  metric = "euclidean",
  min.dist = 0.3,
  spread = 1,
  seed.use = 42L,
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

- graph:

  Optional name of a precomputed Seurat KNN graph produced by
  [`RunWeightedNeighbors()`](https://lachland.github.io/weightedUMAP/reference/RunWeightedNeighbors.md)
  (e.g. `"wt_nn"`). When supplied, UMAP is run on the weighted PC
  embedding that generated the graph (stored by
  [`RunWeightedNeighbors()`](https://lachland.github.io/weightedUMAP/reference/RunWeightedNeighbors.md)
  as `"<prefix>.pca"`, where `<prefix>` is derived by stripping the
  `_nn` suffix from the graph name). This ensures the UMAP topology is
  computed in exactly the same weighted space as the clustering graphs.
  The `dims`, `weight.by`, `reduction`, and `metric` arguments are
  ignored when `graph` is supplied. Set `n.neighbors` to match the
  `k.param` used in
  [`RunWeightedNeighbors()`](https://lachland.github.io/weightedUMAP/reference/RunWeightedNeighbors.md).
  Default: `NULL`.

- reduction.name:

  Name under which to store the new UMAP reduction in `object`. Default:
  `"wt.umap"`.

- reduction.key:

  Column‐name prefix for the UMAP dimensions, e.g. `"wtUMAP_"` produces
  column names `wtUMAP_1`, `wtUMAP_2`, …

- n.neighbors:

  Number of nearest neighbours used by UMAP. Default: `30`.

- n.components:

  Number of UMAP dimensions to compute. Default: `2`.

- metric:

  Distance metric passed to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).
  Default: `"euclidean"`.

- min.dist:

  Effective minimum distance between embedded points. Default: `0.3`.

- spread:

  Effective scale of embedded points. Default: `1`.

- seed.use:

  Random seed for reproducibility. Set to `NULL` to skip
  [`set.seed()`](https://rdrr.io/r/base/Random.html). Default: `42`.

- weight.factor:

  A number in `[0, 1]` that scales the influence of the chosen weighting
  scheme. At `1` (default) the full variance-derived weights are
  applied. At `0` all PCs contribute equally, equivalent to
  `weight.by = "none"` (standard UMAP). Intermediate values blend
  uniformly-weighted and fully-weighted embeddings.

- log.scale:

  Logical. If `TRUE`, [`log1p()`](https://rdrr.io/r/base/Log.html) is
  applied to the weights before the `weight.factor` interpolation,
  compressing the dynamic range so that intermediate PCs receive
  relatively more influence. Ignored when `weight.by = "none"`. Default:
  `FALSE`.

- verbose:

  Print progress messages and per-PC weights? Default: `TRUE`.

- ...:

  Additional arguments forwarded to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).

## Value

The input Seurat `object` with a new `DimReduc` named `reduction.name`
appended. The `misc` slot of this reduction stores:

- `weight.by` — chosen weighting scheme

- `weights` — numeric vector of weights (one per dimension used)

- `source.reduction` — name of the source reduction

- `dims.used` — integer vector of dimensions used

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(wUMAP)
library(patchwork)

# Assumes pbmc has PCA already run

# Standard UMAP — all PCs weighted equally
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "none",
                        reduction.name = "umap.std")

# Weighted UMAP — PCs scaled by proportion of variance explained
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "prop.var",
                        reduction.name = "wt.umap")

# Log-scaled weights — compresses dynamic range so intermediate PCs
# contribute more relative to the dominant PCs
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "prop.var",
                        log.scale = TRUE, reduction.name = "wt.umap.log")

# Standard deviation weighting
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "stdev",
                        reduction.name = "wt.umap.sd")

# Compare standard vs weighted side by side
p1 <- DimPlot(pbmc, reduction = "umap.std",   label = TRUE) +
  ggtitle("Standard UMAP")
p2 <- DimPlot(pbmc, reduction = "wt.umap",    label = TRUE) +
  ggtitle("Weighted UMAP (prop.var)")
p3 <- DimPlot(pbmc, reduction = "wt.umap.log", label = TRUE) +
  ggtitle("Weighted UMAP (prop.var, log.scale)")
p1 | p2
} # }
```
