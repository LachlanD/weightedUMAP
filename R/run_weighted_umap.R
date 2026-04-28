#' Run Variance-Weighted UMAP on a Seurat Object
#'
#' Scales each principal component embedding by a weight derived from the
#' variance explained by that component, then runs UMAP via \pkg{uwot}.
#' The result is stored as a new `DimReduc` in the returned Seurat object.
#'
#' @param object A Seurat object with a PCA (or other linear) reduction
#'   already computed.
#' @param reduction Name of the dimensional reduction to use as input.
#'   Default: `"pca"`.
#' @param dims Integer vector of dimensions to use. If `NULL` (default), all
#'   available dimensions are used and a warning is raised when there are more
#'   than 50.
#' @param weight.by Weighting scheme applied to PC scores before UMAP.
#'   One of:
#'   \describe{
#'     \item{`"stdev"`}{Standard deviation, normalised (`sdev / sum(sdev)`).
#'       Default. Gently up-weights early PCs while keeping intermediate ones
#'       in play, giving a good balance between signal emphasis and layout
#'       completeness.}
#'     \item{`"prop.var"`}{Proportion of variance explained
#'       (`sdev^2 / sum(sdev^2)`); more aggressively up-weights the dominant
#'       PCs.  Use with `weight.factor < 1` or `log.scale = TRUE` to avoid
#'       over-compression on datasets where PC 1 explains a large fraction of
#'       variance.}
#'     \item{`"none"`}{No weighting; equivalent to standard UMAP on PCA
#'       scores.}
#'   }
#' @param mp.filter Logical. If `TRUE`, PCs whose variance is at or below the
#'   Marchenko-Pastur bulk-noise threshold
#'   (\eqn{\lambda_{\max} = (1 + \sqrt{p/n})^2}) are zeroed out after the
#'   chosen `weight.by` weights are computed, and the remaining weights are
#'   renormalised.  This makes MP an orthogonal filtering step compatible with
#'   any weighting scheme.  Ignored when `weight.by = "none"` and no PCs are
#'   above the threshold.  Default: `FALSE`.
#' @param reduction.name Name under which to store the new UMAP reduction in
#'   `object`. Default: `"wt.umap"`.
#' @param reduction.key Column‐name prefix for the UMAP dimensions, e.g.
#'   `"wtUMAP_"` produces column names `wtUMAP_1`, `wtUMAP_2`, …
#' @param n.neighbors Number of nearest neighbours used by UMAP. Default: `30`.
#' @param n.components Number of UMAP dimensions to compute. Default: `2`.
#' @param metric Distance metric passed to [uwot::umap()]. Default:
#'   `"euclidean"`.
#' @param min.dist Effective minimum distance between embedded points.
#'   Default: `0.3`.
#' @param spread Effective scale of embedded points. Default: `1`.
#' @param seed.use Random seed for reproducibility. Set to `NULL` to skip
#'   `set.seed()`. Default: `42`.
#' @param weight.factor A number in `[0, 1]` that scales the influence of the
#'   chosen weighting scheme.  At `1` (default) the full variance-derived
#'   weights are applied.  At `0` all PCs contribute equally, equivalent to
#'   `weight.by = "none"` (standard UMAP).  Intermediate values blend
#'   uniformly-weighted and fully-weighted embeddings.
#' @param log.scale Logical. If `TRUE`, `log1p()` is applied to the weights
#'   before the `weight.factor` interpolation, compressing the dynamic range
#'   so that intermediate PCs receive relatively more influence.  Ignored
#'   when `weight.by = "none"`.  Default: `FALSE`.
#' @param verbose Print progress messages and per-PC weights? Default: `TRUE`.
#' @param graph Optional name of a precomputed Seurat KNN graph produced by
#'   [RunWeightedNeighbors()] (e.g. `"wt_nn"`).  When supplied, UMAP is run
#'   on the weighted PC embedding that generated the graph (stored by
#'   `RunWeightedNeighbors()` as `"<prefix>.pca"`, where `<prefix>` is derived
#'   by stripping the `_nn` suffix from the graph name).  This ensures the
#'   UMAP topology is computed in exactly the same weighted space as the
#'   clustering graphs.  The `dims`, `weight.by`, `reduction`, and `metric`
#'   arguments are ignored when `graph` is supplied.  Set `n.neighbors` to
#'   match the `k.param` used in [RunWeightedNeighbors()].  Default: `NULL`.
#' @param ... Additional arguments forwarded to [uwot::umap()].
#'
#' @return The input Seurat `object` with a new `DimReduc` named
#'   `reduction.name` appended. The `misc` slot of this reduction stores:
#'   * `weight.by` — chosen weighting scheme
#'   * `weights`    — numeric vector of weights (one per dimension used)
#'   * `source.reduction` — name of the source reduction
#'   * `dims.used` — integer vector of dimensions used
#'
#' @importFrom SeuratObject Embeddings Stdev DefaultAssay CreateDimReducObject
#' @importFrom uwot umap
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(wUMAP)
#' library(patchwork)
#'
#' # Assumes pbmc has PCA already run
#'
#' # Weighted UMAP — PCs scaled by standard deviation (default, recommended)
#' pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, reduction.name = "wt.umap")
#'
#' # Standard UMAP for comparison — all PCs weighted equally
#' pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "none",
#'                         reduction.name = "umap.std")
#'
#' # MP filtering: zero out noise PCs, then apply stdev weights to signal PCs
#' pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, mp.filter = TRUE,
#'                         reduction.name = "wt.umap.mp")
#'
#' # Compare standard vs weighted
#' p1 <- DimPlot(pbmc, reduction = "umap.std", label = TRUE) +
#'   ggtitle("Standard UMAP")
#' p2 <- DimPlot(pbmc, reduction = "wt.umap",  label = TRUE) +
#'   ggtitle("Weighted UMAP (stdev)")
#' p1 | p2
#' }
RunWeightedUMAP <- function(
    object,
    reduction      = "pca",
    dims           = NULL,
    weight.by      = c("stdev", "prop.var", "none"),
    graph          = NULL,
    reduction.name = "wt.umap",
    reduction.key  = "wtUMAP_",
    n.neighbors    = 30L,
    n.components   = 2L,
    metric         = "euclidean",
    min.dist       = 0.3,
    spread         = 1,
    seed.use       = 42L,
    weight.factor  = 1,
    log.scale      = FALSE,
    mp.filter      = FALSE,
    verbose        = TRUE,
    ...
) {
  # ── Input validation ────────────────────────────────────────────────────────
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop(
      "The 'uwot' package is required. Install it with:\n",
      "  install.packages('uwot')",
      call. = FALSE
    )
  }

  weight.by <- match.arg(weight.by)

  if (!is.numeric(weight.factor) || length(weight.factor) != 1 ||
      weight.factor < 0 || weight.factor > 1) {
    stop("'weight.factor' must be a single number between 0 and 1.",
         call. = FALSE)
  }

  if (!is.logical(log.scale) || length(log.scale) != 1) {
    stop("'log.scale' must be TRUE or FALSE.", call. = FALSE)
  }

  # ── Graph-based path (uses weighted embedding that produced the graph) ──────
  if (!is.null(graph)) {
    if (!graph %in% names(object@graphs)) {
      stop(
        sprintf(
          "Graph '%s' not found in object. Available graphs: %s",
          graph,
          paste(names(object@graphs), collapse = ", ")
        ),
        call. = FALSE
      )
    }

    # Derive source weighted reduction by convention: "wt_nn" → "wt.pca"
    if (grepl("_snn$", graph)) {
      stop(
        sprintf(
          "Graph '%s' is an SNN graph and cannot be used directly for UMAP. ",
          graph
        ),
        "Pass the corresponding KNN graph instead (replace '_snn' with '_nn').",
        call. = FALSE
      )
    }
    prefix          <- sub("_nn$", "", graph)
    wt_red_name     <- paste0(prefix, ".pca")

    if (!wt_red_name %in% names(object@reductions)) {
      stop(
        sprintf(
          paste0(
            "Could not find source weighted reduction '%s' (expected to be ",
            "created by RunWeightedNeighbors() alongside graph '%s'). ",
            "Available reductions: %s"
          ),
          wt_red_name, graph,
          paste(names(object@reductions), collapse = ", ")
        ),
        call. = FALSE
      )
    }

    wt_emb     <- Embeddings(object[[wt_red_name]])
    cell_names <- rownames(wt_emb)

    if (verbose) {
      message(sprintf(
        "[wUMAP] Using weighted embedding '%s' (%d cells x %d dims, k = %d)",
        wt_red_name, nrow(wt_emb), ncol(wt_emb), n.neighbors
      ))
    }

    if (!is.null(seed.use)) set.seed(seed.use)

    umap_mat <- uwot::umap(
      X            = wt_emb,
      n_neighbors  = as.integer(n.neighbors),
      n_components = as.integer(n.components),
      metric       = metric,
      min_dist     = min.dist,
      spread       = spread,
      verbose      = verbose,
      ...
    )

    rownames(umap_mat) <- cell_names
    colnames(umap_mat) <- paste0(reduction.key, seq_len(n.components))

    source_assay <- tryCatch(
      DefaultAssay(object[[wt_red_name]]),
      error = function(e) DefaultAssay(object)
    )

    object[[reduction.name]] <- CreateDimReducObject(
      embeddings = umap_mat,
      key        = reduction.key,
      assay      = source_assay,
      misc       = list(source.graph = graph)
    )

    if (verbose) {
      message(sprintf(
        "[wUMAP] Stored as reduction '%s'. Access via Embeddings(object[['%s']]).",
        reduction.name, reduction.name
      ))
    }

    return(object)
  }

  if (!reduction %in% names(object@reductions)) {
    stop(
      sprintf(
        "Reduction '%s' not found in object. Available reductions: %s",
        reduction,
        paste(names(object@reductions), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # ── Compute weighted embeddings via shared helper ────────────────────────────
  wt           <- .compute_weighted_embeddings(object, reduction, dims,
                                               weight.by, weight.factor,
                                               log.scale, verbose,
                                               mp.filter = mp.filter)
  emb          <- wt$emb
  weighted_emb <- wt$weighted_emb
  dims         <- wt$dims

  # ── Run UMAP ─────────────────────────────────────────────────────────────────
  if (!is.null(seed.use)) set.seed(seed.use)

  umap_mat <- uwot::umap(
    X           = weighted_emb,
    n_neighbors = as.integer(n.neighbors),
    n_components = as.integer(n.components),
    metric      = metric,
    min_dist    = min.dist,
    spread      = spread,
    verbose     = verbose,
    ...
  )

  # ── Name rows / columns ──────────────────────────────────────────────────────
  rownames(umap_mat) <- rownames(emb)
  colnames(umap_mat) <- paste0(reduction.key, seq_len(n.components))

  # ── Build DimReduc and store ─────────────────────────────────────────────────
  source_assay <- tryCatch(
    DefaultAssay(object[[reduction]]),
    error = function(e) DefaultAssay(object)
  )

  dimreduc <- CreateDimReducObject(
    embeddings = umap_mat,
    key        = reduction.key,
    assay      = source_assay,
    misc       = list(
      weight.by        = weight.by,
      mp.filter        = mp.filter,
      weights          = wt$weights,
      source.reduction = reduction,
      dims.used        = dims
    )
  )

  object[[reduction.name]] <- dimreduc

  if (verbose) {
    message(sprintf(
      "[wUMAP] Stored as reduction '%s'. Access via Embeddings(object[['%s']]).",
      reduction.name, reduction.name
    ))
  }

  return(object)
}
