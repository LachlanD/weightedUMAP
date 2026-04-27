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
#'     \item{`"pct.var"`}{Percentage of variance explained
#'       (`sdev^2 / sum(sdev^2) * 100`). Default.}
#'     \item{`"prop.var"`}{Proportion of variance explained
#'       (`sdev^2 / sum(sdev^2)`).}
#'     \item{`"eigenvalue"`}{Eigenvalue (`sdev^2`).}
#'     \item{`"stdev"`}{Standard deviation (`sdev`).}
#'     \item{`"none"`}{No weighting; equivalent to standard UMAP on PCA
#'       scores.}
#'   }
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
#' @param verbose Print progress messages and per-PC weights? Default: `TRUE`.
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
#'
#' # Assumes pbmc has PCA already run
#' pbmc <- RunWeightedUMAP(pbmc, dims = 1:30, weight.by = "pct.var")
#' DimPlot(pbmc, reduction = "wt.umap")
#' }
RunWeightedUMAP <- function(
    object,
    reduction      = "pca",
    dims           = NULL,
    weight.by      = c("pct.var", "prop.var", "eigenvalue", "stdev", "none"),
    reduction.name = "wt.umap",
    reduction.key  = "wtUMAP_",
    n.neighbors    = 30L,
    n.components   = 2L,
    metric         = "euclidean",
    min.dist       = 0.3,
    spread         = 1,
    seed.use       = 42L,
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

  # ── Extract embeddings & standard deviations ────────────────────────────────
  emb_full  <- Embeddings(object[[reduction]])   # cells × dims
  sdev_full <- Stdev(object[[reduction]])         # numeric vector, length = dims

  if (length(sdev_full) == 0) {
    stop(
      sprintf(
        "No standard deviations found for reduction '%s'. ",
        reduction
      ),
      "Ensure you ran PCA with Seurat (e.g. RunPCA()) so that stdev values ",
      "are stored.",
      call. = FALSE
    )
  }

  n_dims_available <- min(ncol(emb_full), length(sdev_full))

  if (is.null(dims)) {
    if (n_dims_available > 50) {
      warning(
        sprintf(
          "'dims' not specified and %d dimensions are available. ",
          n_dims_available
        ),
        "Using all of them. Consider supplying 'dims' explicitly, e.g. ",
        "dims = 1:30.",
        call. = FALSE
      )
    }
    dims <- seq_len(n_dims_available)
  } else {
    dims <- as.integer(dims)
    bad  <- dims[dims < 1 | dims > n_dims_available]
    if (length(bad) > 0) {
      stop(
        sprintf(
          "Requested dims (%s) exceed available range (1:%d).",
          paste(bad, collapse = ", "),
          n_dims_available
        ),
        call. = FALSE
      )
    }
  }

  emb  <- emb_full[, dims, drop = FALSE]
  sdev <- sdev_full[dims]

  # ── Compute weights ──────────────────────────────────────────────────────────
  weights <- switch(
    weight.by,
    pct.var    = sdev^2 / sum(sdev^2) * 100,
    prop.var   = sdev^2 / sum(sdev^2),
    eigenvalue = sdev^2,
    stdev      = sdev,
    none       = rep(1.0, length(sdev))
  )

  if (verbose) {
    message(sprintf("[wUMAP] Weight scheme : %s", weight.by))
    message(sprintf("[wUMAP] Dimensions    : %d (dims %d–%d)",
                    length(dims), min(dims), max(dims)))
    max_show <- min(length(weights), 10L)
    wt_str   <- paste(
      sprintf("PC%d=%.3g", dims[seq_len(max_show)], weights[seq_len(max_show)]),
      collapse = ", "
    )
    if (length(weights) > 10L) wt_str <- paste0(wt_str, ", ...")
    message(sprintf("[wUMAP] Weights       : %s", wt_str))
  }

  # ── Scale embeddings ─────────────────────────────────────────────────────────
  weighted_emb <- sweep(emb, 2L, weights, `*`)

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
      weights          = weights,
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
