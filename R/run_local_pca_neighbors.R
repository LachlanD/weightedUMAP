#' Build a Local-PCA Nearest-Neighbour Graph
#'
#' Computes per-cell local PCA distances (see \code{\link{RunLocalPCAUMAP}})
#' and stores the resulting KNN and SNN graphs in the Seurat object.  The
#' graphs can be passed to \code{\link[Seurat]{FindClusters}} (via the
#' \code{graph.name} argument) and to \code{\link{RunLocalPCAUMAP}} (via the
#' \code{nn_method} argument in \code{...}), ensuring that clustering and
#' visualisation share exactly the same local-PCA nearest-neighbour structure.
#'
#' The SNN graph is a sparse matrix of Jaccard similarities, pruned at
#' \code{prune.SNN} (default \eqn{1/15}), matching Seurat's convention.
#'
#' @param object A Seurat object with PCA already run.
#' @param reduction Name of the PCA reduction to use. Default: \code{"pca"}.
#' @param dims Integer vector of PC dimensions to use. Default: all available.
#' @param k.param Number of nearest neighbours. Default: \code{20L}.
#' @param local.dims Number of local PCA directions to retain when computing
#'   distances.  Must be \eqn{\leq k.param}.
#'   Default: \code{NULL} (uses \code{k.param} directions).
#' @param local.weight.by Weighting scheme for local PC directions.  See
#'   \code{\link{RunLocalPCAUMAP}} for full descriptions.  One of
#'   \code{"stdev"} (default), \code{"prop.var"}, \code{"none"}.
#' @param mp.filter Logical. If \code{TRUE}, local PC directions at or below
#'   the per-neighbourhood Marchenko-Pastur noise threshold are zeroed out
#'   after the base weights are computed.  Can be combined with any
#'   \code{local.weight.by} scheme.  Default: \code{FALSE}.
#' @param prune.SNN Jaccard similarity threshold below which edges are removed
#'   from the SNN graph. Default: \code{1/15}.
#' @param prefix Prefix for the two stored graphs: \code{<prefix>_nn} (KNN,
#'   binary) and \code{<prefix>_snn} (SNN, Jaccard-weighted).  Pass
#'   \code{<prefix>_snn} to \code{FindClusters(graph.name = ...)}.
#'   Default: \code{"lp"}.
#' @param verbose Print progress messages. Default: \code{TRUE}.
#'
#' @return The input Seurat \code{object} with two new Graph objects:
#'   \itemize{
#'     \item \code{<prefix>_nn} — binary KNN graph (symmetric).
#'     \item \code{<prefix>_snn} — Jaccard-weighted SNN graph (sparse,
#'           pruned).
#'   }
#'
#' @seealso \code{\link{RunLocalPCAUMAP}}, \code{\link{RunWeightedNeighbors}}
#'
#' @importFrom SeuratObject Embeddings Stdev DefaultAssay as.Graph
#' @importFrom Matrix sparseMatrix tcrossprod drop0 t
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(wUMAP)
#'
#' # Assumes pbmc has PCA already run
#'
#' # 1. Build local PCA KNN/SNN graphs — stdev weighting (default, recommended)
#' pbmc <- RunLocalPCANeighbors(pbmc, dims = 1:30, k.param = 20,
#'                              prefix = "lp")
#'
#' # 2. Cluster on the local PCA SNN graph
#' pbmc <- FindClusters(pbmc, graph.name = "lp_snn", resolution = 0.5)
#'
#' # 3. Visualise with matching local PCA UMAP (same k, same weighting)
#' pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 20,
#'                         reduction.name = "lp.umap")
#'
#' DimPlot(pbmc, reduction = "lp.umap", label = TRUE)
#'
#' # With MP filtering: discard local noise directions, then weight by stdev
#' pbmc <- RunLocalPCANeighbors(pbmc, dims = 1:30, k.param = 20,
#'                              mp.filter = TRUE, prefix = "lp.mp")
#' pbmc <- FindClusters(pbmc, graph.name = "lp.mp_snn", resolution = 0.5)
#' pbmc <- RunLocalPCAUMAP(pbmc, dims = 1:30, k.param = 20,
#'                         mp.filter = TRUE, reduction.name = "lp.umap.mp")
#' }
RunLocalPCANeighbors <- function(
    object,
    reduction       = "pca",
    dims            = NULL,
    k.param         = 20L,
    local.dims      = NULL,
    local.weight.by = c("stdev", "prop.var", "none"),
    mp.filter       = FALSE,
    prune.SNN       = 1 / 15,
    prefix          = "lp",
    verbose         = TRUE
) {
  if (!requireNamespace("RANN", quietly = TRUE))
    stop("Package 'RANN' is required. Install with install.packages('RANN').",
         call. = FALSE)

  local.weight.by <- match.arg(local.weight.by)

  if (!is.numeric(prune.SNN) || length(prune.SNN) != 1 ||
      prune.SNN < 0 || prune.SNN > 1)
    stop("'prune.SNN' must be a single number in [0, 1].", call. = FALSE)

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
      "[wUMAP] RunLocalPCANeighbors: %d cells x %d dims, k = %d, local.dims = %d, local.weight.by = '%s', mp.filter = %s",
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
  nn_idx     <- res$nn_idx     # n × k (self excluded)
  local_dist <- res$local_dist # n × k

  cell_names <- rownames(emb)

  # ── Build KNN graph (binary, symmetric) ───────────────────────────────────
  if (verbose) message("[wUMAP] Building KNN graph ...")

  i_knn_raw <- rep(seq_len(n), each = k.param)
  j_knn_raw <- as.vector(t(nn_idx))

  # Drop any entries where RANN placed a cell as its own k-NN (beyond column 1)
  not_self  <- i_knn_raw != j_knn_raw
  i_knn     <- i_knn_raw[not_self]
  j_knn     <- j_knn_raw[not_self]

  knn_dir <- Matrix::sparseMatrix(
    i = i_knn, j = j_knn, x = 1L,
    dims = c(n, n),
    dimnames = list(cell_names, cell_names)
  )
  # Symmetrize: union of directed edges → binary
  knn_sym      <- knn_dir + Matrix::t(knn_dir)
  knn_sym@x[]  <- 1L

  # ── Build SNN graph (Jaccard-weighted) ────────────────────────────────────
  if (verbose) message("[wUMAP] Building SNN graph (Jaccard) ...")

  # Binary adjacency WITH self for shared-neighbor counting.
  # Build off-diagonal part from (deduplicated) knn edges, then add Diagonal(n)
  # so each cell's self-entry is exactly 1, regardless of RANN quirks.
  k1 <- k.param + 1L
  A   <- Matrix::sparseMatrix(i = i_knn, j = j_knn, x = 1L, dims = c(n, n)) +
         Matrix::Diagonal(n)

  # Shared-neighbor count matrix
  shared    <- Matrix::tcrossprod(A)             # n×n, [i,j] = |N(i) ∩ N(j)|
  shared_ij <- shared[cbind(i_knn, j_knn)]
  jaccard   <- shared_ij / (2L * k1 - shared_ij)

  # Symmetric SNN via unique unordered pairs.
  # Jaccard is symmetric: shared[i,j] = shared[j,i] for all pairs.
  # For mutual kNN edges (i→j AND j→i) both yield the same Jaccard value;
  # deduplicate to avoid storing 2× before building the symmetric matrix.
  pair_lo  <- pmin(i_knn, j_knn)
  pair_hi  <- pmax(i_knn, j_knn)
  uniq_idx <- !duplicated(cbind(pair_lo, pair_hi))

  ui   <- pair_lo[uniq_idx]
  uj   <- pair_hi[uniq_idx]
  ujac <- jaccard[uniq_idx]

  snn_sym <- Matrix::sparseMatrix(
    i = c(ui, uj),
    j = c(uj, ui),
    x = c(ujac, ujac),
    dims = c(n, n),
    dimnames = list(cell_names, cell_names)
  )

  diag(snn_sym)      <- 0
  snn_sym            <- Matrix::drop0(snn_sym)
  snn_sym[snn_sym < prune.SNN] <- 0
  snn_sym            <- Matrix::drop0(snn_sym)

  # ── Store as Seurat Graph objects ─────────────────────────────────────────
  source_assay <- tryCatch(
    DefaultAssay(object[[reduction]]),
    error = function(e) DefaultAssay(object)
  )

  knn_graph <- SeuratObject::as.Graph(knn_sym)
  snn_graph <- SeuratObject::as.Graph(snn_sym)

  tryCatch({
    SeuratObject::DefaultAssay(knn_graph) <- source_assay
    SeuratObject::DefaultAssay(snn_graph) <- source_assay
  }, error = function(e) NULL)

  object[[paste0(prefix, "_nn")]]  <- knn_graph
  object[[paste0(prefix, "_snn")]] <- snn_graph

  if (verbose)
    message(sprintf(
      "[wUMAP] Stored graphs '%s_nn' and '%s_snn'.", prefix, prefix
    ))

  object
}
