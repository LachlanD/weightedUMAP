# Run PC-Weighted UMAP on a Seurat Object

Runs UMAP via uwot with an optionally modified PC-space distance: each
principal component axis can be scaled by a variance-derived weight
(heuristic) and/or noise PCs can be zeroed by the Marchenko-Pastur
criterion (`mp.filter`). The result is stored as a new `DimReduc` in the
returned Seurat object.

## Usage

``` r
RunWeightedUMAP(
  object,
  reduction = "pca",
  dims = NULL,
  weight.by = c("stdev", "prop.var", "none"),
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
  mp.filter = FALSE,
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

  `"stdev"`

  :   Standard deviation, normalised (`sdev / sum(sdev)`). Default.
      Mildest transformation — early PCs receive somewhat more weight.
      Whether this improves results is dataset-dependent.

  `"prop.var"`

  :   Proportion of variance explained (`sdev^2 / sum(sdev^2)`); more
      aggressively up-weights the dominant PCs. Use with
      `weight.factor < 1` or `log.scale = TRUE` to avoid
      over-compression on datasets where PC 1 explains a large fraction
      of variance.

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

- mp.filter:

  Logical. If `TRUE`, PCs whose variance is at or below the
  Marchenko-Pastur bulk-noise threshold (\\\lambda\_{\max} = (1 +
  \sqrt{p/n})^2\\) are zeroed out after the chosen `weight.by` weights
  are computed, and the remaining weights are renormalised. This makes
  MP an orthogonal filtering step compatible with any weighting scheme.
  Ignored when `weight.by = "none"` and no PCs are above the threshold.
  Default: `FALSE`.

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

# Weighted UMAP — PCs scaled by standard deviation (default, recommended)
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, reduction.name = "wt.umap")

# Standard UMAP for comparison — all PCs weighted equally
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "none",
                        reduction.name = "umap.std")

# MP filtering: zero out noise PCs, then apply stdev weights to signal PCs
pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, mp.filter = TRUE,
                        reduction.name = "wt.umap.mp")

# Compare standard vs weighted
p1 <- DimPlot(pbmc, reduction = "umap.std", label = TRUE) +
  ggtitle("Standard UMAP")
p2 <- DimPlot(pbmc, reduction = "wt.umap",  label = TRUE) +
  ggtitle("Weighted UMAP (stdev)")
p1 | p2
} # }
```
