#' UMAP with Local PCA Distances
#'
#' Computes UMAP using a **local PCA metric**: for each cell, the distance to
#' each of its \eqn{k} nearest neighbours is re-measured in a locally-fitted
#' PCA basis built from that neighbourhood.  This captures anisotropic local
#' structure (e.g. the tangent direction of a trajectory) and de-emphasises
#' transverse noise, producing a richer distance measure than a single global
#' Euclidean metric.
#'
#' Optionally, each local PC direction can be weighted by its local variance
#' contribution (\code{local.weight.by}), analogous to how \code{weight.by}
#' weights global PCs in \code{\link{RunWeightedUMAP}}.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Find the global \code{k.param} nearest neighbours of every cell in
#'     the supplied PCA embedding (using \pkg{RANN}).
#'   \item For each cell \eqn{i}, collect its \eqn{k}-neighbourhood, centre by
#'     the neighbourhood mean, and compute a compact SVD.
#'   \item Optionally weight each of the \code{local.dims} local PC directions
#'     by its local variance contribution (\code{local.weight.by}).
#'   \item Report the (weighted) Euclidean norm in local-PC space as the
#'     refined distance.
#'   \item Pass the \eqn{n \times k} index/distance matrices directly to
#'     \code{\link[uwot]{umap}} as a precomputed k-NN graph.
#' }
#'
#' @param object A Seurat object with PCA already run.
#' @param reduction Name of the PCA reduction to use. Default: \code{"pca"}.
#' @param dims Integer vector of PC dimensions to use. Default: all available.
#' @param k.param Number of nearest neighbours for the local PCA.
#'   Default: \code{30L}.
#' @param local.dims Number of local PCA directions to retain when computing
#'   refined distances.  Must be \eqn{\leq} \code{k.param}.
#'   Default: \code{NULL} (uses \code{k.param} directions).
#' @param local.weight.by Weighting scheme applied to local PC directions.
#'   Analogous to \code{weight.by} in \code{\link{RunWeightedUMAP}}, but
#'   applied to the \emph{local} singular values of each neighbourhood SVD.
#'   \describe{
#'     \item{\code{"stdev"}}{Weight by local standard deviation
#'       (\eqn{s_l / \sum s}).  Default.  Gently emphasises the dominant
#'       local direction without collapsing minor variation.}
#'     \item{\code{"prop.var"}}{Weight each local PC by its proportion of
#'       local variance (\eqn{s_l^2 / \sum s^2}); more aggressively
#'       emphasises the leading local direction.}
#'     \item{\code{"none"}}{No weighting; equal contribution from all local
#'       PCs.}
#'   }
#' @param mp.filter Logical. If \code{TRUE}, local PC directions whose local
#'   variance is at or below the per-neighbourhood Marchenko-Pastur noise
#'   threshold (\eqn{\lambda_{\max} = (1 + \sqrt{d/k})^2}) are zeroed out
#'   after the \code{local.weight.by} weights are computed, and the remaining
#'   weights are renormalised.  Can be combined with any \code{local.weight.by}
#'   scheme.  Default: \code{FALSE}.
#' @param reduction.name Name under which to store the UMAP result.
#'   Default: \code{"lp.umap"}.
#' @param reduction.key Column-name prefix for the UMAP dimensions.
#'   Default: \code{"lpUMAP_"}.
#' @param n.components Number of UMAP dimensions. Default: \code{2L}.
#' @param min.dist UMAP \code{min_dist} parameter. Default: \code{0.3}.
#' @param spread UMAP \code{spread} parameter. Default: \code{1}.
#' @param seed.use Random seed. Default: \code{42L}.
#' @param verbose Print progress messages. Default: \code{TRUE}.
#' @param ... Additional arguments forwarded to \code{\link[uwot]{umap}}.
#'
#' @return The input Seurat \code{object} with a new DimReduc stored under
#'   \code{reduction.name}.  The \code{misc} slot records
#'   \code{source.reduction}, \code{dims.used}, \code{k.param},
#'   \code{local.dims}, and \code{local.weight.by}.
#'
#' @seealso \code{\link{RunLocalPCANeighbors}}, \code{\link{RunWeightedUMAP}},
#'   \code{\link{RunWeightedNeighbors}}
#'
#' @importFrom SeuratObject Embeddings Stdev DefaultAssay CreateDimReducObject
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(wUMAP)
#'
#' # Assumes pbmc has PCA already run
#'
#' # Local PCA UMAP with stdev weighting of local PC directions (default)
#' pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 30,
#'                         reduction.name = "lp.umap")
#'
#' # Unweighted local PCA for comparison
#' pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 30,
#'                         local.weight.by = "none",
#'                         reduction.name = "lp.umap.unweighted")
#'
#' # With MP filtering: zero local noise directions, then apply stdev weights
#' pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 30,
#'                         mp.filter = TRUE,
#'                         reduction.name = "lp.umap.mp")
#'
#' # Compare standard, unweighted local PCA, and weighted local PCA
#' p1 <- DimPlot(pbmc, reduction = "umap",               label = TRUE) +
#'   ggtitle("Standard")
#' p2 <- DimPlot(pbmc, reduction = "lp.umap.unweighted", label = TRUE) +
#'   ggtitle("Local PCA (unweighted)")
#' p3 <- DimPlot(pbmc, reduction = "lp.umap",            label = TRUE) +
#'   ggtitle("Local PCA (stdev)")
#' p1 | p2 | p3
#' }
RunLocalPCAUMAP <- function(
    object,
    reduction       = "pca",
    dims            = NULL,
    k.param         = 30L,
    local.dims      = NULL,
    local.weight.by = c("stdev", "prop.var", "none"),
    mp.filter       = FALSE,
    reduction.name  = "lp.umap",
    reduction.key   = "lpUMAP_",
    n.components    = 2L,
    min.dist        = 0.3,
    spread          = 1,
    seed.use        = 42L,
    verbose         = TRUE,
    ...
) {
  if (!requireNamespace("uwot", quietly = TRUE))
    stop("Package 'uwot' is required. Install with install.packages('uwot').",
         call. = FALSE)
  if (!requireNamespace("RANN", quietly = TRUE))
    stop("Package 'RANN' is required. Install with install.packages('RANN').",
         call. = FALSE)

  local.weight.by <- match.arg(local.weight.by)

  if (!reduction %in% names(object@reductions))
    stop(sprintf("Reduction '%s' not found. Available: %s",
                 reduction, paste(names(object@reductions), collapse = ", ")),
         call. = FALSE)

  emb_full <- Embeddings(object[[reduction]])
  sdev_full <- Stdev(object[[reduction]])
  n_avail  <- min(ncol(emb_full), length(sdev_full))

  if (is.null(dims)) {
    dims <- seq_len(n_avail)
  } else {
    dims <- as.integer(dims)
    bad  <- dims[dims < 1 | dims > n_avail]
    if (length(bad))
      stop(sprintf("Requested dims (%s) exceed available range (1:%d).",
                   paste(bad, collapse = ", "), n_avail), call. = FALSE)
  }

  emb <- emb_full[, dims, drop = FALSE]
  n   <- nrow(emb)
  d   <- ncol(emb)

  k.param <- as.integer(k.param)
  if (k.param < 2L || k.param >= n)
    stop("'k.param' must be between 2 and n-1.", call. = FALSE)

  if (is.null(local.dims)) {
    local.dims <- k.param
  } else {
    local.dims <- as.integer(local.dims)
    if (local.dims < 1L || local.dims > k.param)
      stop("'local.dims' must be between 1 and k.param.", call. = FALSE)
  }
  local.dims <- min(local.dims, d)

  if (verbose)
    message(sprintf(
      "[wUMAP] RunLocalPCAUMAP: %d cells x %d dims, k = %d, local.dims = %d, local.weight.by = '%s', mp.filter = %s",
      n, d, k.param, local.dims, local.weight.by, mp.filter
    ))

  # ── Compute local PCA distances (shared helper) ────────────────────────────
  res <- .compute_local_pca_distances(
    emb             = emb,
    k.param         = k.param,
    local.dims      = local.dims,
    local.weight.by = local.weight.by,
    mp.filter       = mp.filter,
    verbose         = verbose
  )
  nn_idx     <- res$nn_idx
  local_dist <- res$local_dist

  if (verbose) message("[wUMAP] Step 3/3: running UMAP on local PCA graph ...")

  # ── UMAP with precomputed kNN ─────────────────────────────────────────────
  if (!is.null(seed.use)) set.seed(seed.use)

  umap_mat <- uwot::umap(
    X            = NULL,
    nn_method    = list(idx = nn_idx, dist = local_dist),
    n_neighbors  = k.param,
    n_components = as.integer(n.components),
    min_dist     = min.dist,
    spread       = spread,
    verbose      = verbose,
    ...
  )

  rownames(umap_mat) <- rownames(emb)
  colnames(umap_mat) <- paste0(reduction.key, seq_len(n.components))

  source_assay <- tryCatch(
    DefaultAssay(object[[reduction]]),
    error = function(e) DefaultAssay(object)
  )

  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = umap_mat,
    key        = reduction.key,
    assay      = source_assay,
    misc       = list(
      source.reduction = reduction,
      dims.used        = dims,
      k.param          = k.param,
      local.dims       = local.dims,
      local.weight.by  = local.weight.by,
      mp.filter        = mp.filter
    )
  )

  if (verbose)
    message(sprintf("[wUMAP] Stored as reduction '%s'.", reduction.name))

  object
}
