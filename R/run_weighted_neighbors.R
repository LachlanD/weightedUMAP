#' Build a Variance-Weighted Nearest-Neighbour Graph
#'
#' Scales each principal component by a variance-derived weight (the same
#' weighting as \code{\link{RunWeightedUMAP}}), stores the result as a new
#' \code{DimReduc}, and builds a KNN/SNN graph via Seurat's
#' \code{FindNeighbors()}.  The resulting graphs can be passed directly to
#' \code{\link[Seurat]{FindClusters}} and, via the \code{graph} argument of
#' \code{\link{RunWeightedUMAP}}, for UMAP — ensuring that clustering and
#' visualisation share exactly the same nearest-neighbour structure.
#'
#' @inheritParams RunWeightedUMAP
#' @param k.param Number of nearest neighbours for the KNN graph.
#'   Default: \code{20L}.
#' @param reduction.name Name under which to store the intermediate weighted-PC
#'   embedding in \code{object}.  Default: \code{"wt.pca"}.
#' @param reduction.key Column-name prefix for the weighted-PC dimensions.
#'   Default: \code{"wtPCA_"}.
#' @param graph.name Prefix for the two graphs stored in the Seurat object:
#'   \code{<graph.name>_nn} (KNN, used for UMAP) and
#'   \code{<graph.name>_snn} (SNN, used for clustering).
#'   Pass this prefix to \code{FindClusters(graph.name = "<graph.name>_snn")}
#'   and \code{RunWeightedUMAP(graph = "<graph.name>_nn")}.
#'   Default: \code{"wt"}.
#' @param weight.factor A number in \code{[0, 1]} controlling the influence of
#'   the weighting scheme.  At \code{1} (default) full variance-derived weights
#'   are applied.  At \code{0} all PCs contribute equally (standard
#'   nearest-neighbour graph).  Must match the \code{weight.factor} used in the
#'   subsequent \code{\link{RunWeightedUMAP}} call to keep clustering and
#'   visualisation in the same space.
#' @param ... Additional arguments forwarded to
#'   \code{\link[Seurat]{FindNeighbors}}.
#'
#' @return The input Seurat \code{object} with:
#'   \itemize{
#'     \item A new \code{DimReduc} stored under \code{reduction.name}
#'           containing the variance-weighted PC embeddings.
#'     \item Two neighbour graphs (\code{<graph.name>_nn} and
#'           \code{<graph.name>_snn}) suitable for downstream
#'           \code{FindClusters()} and \code{RunWeightedUMAP(graph = ...)}.
#'   }
#'
#' @seealso \code{\link{RunWeightedUMAP}}
#'
#' @importFrom SeuratObject Embeddings Stdev DefaultAssay CreateDimReducObject
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(wUMAP)
#'
#' # 1. Compute the weighted KNN/SNN graphs (same space for clustering + UMAP)
#' pbmc <- RunWeightedNeighbors(pbmc, dims = 1:30, weight.by = "prop.var",
#'                              graph.name = "wt")
#'
#' # 2. Cluster on the weighted SNN graph
#' pbmc <- FindClusters(pbmc, graph.name = "wt_snn")
#'
#' # 3. UMAP from the same weighted KNN graph — topology is identical to step 2
#' pbmc <- RunWeightedUMAP(pbmc, graph = "wt_nn", n.neighbors = 20,
#'                         reduction.name = "wt.umap")
#'
#' DimPlot(pbmc, reduction = "wt.umap", label = TRUE)
#' }
RunWeightedNeighbors <- function(
    object,
    reduction      = "pca",
    dims           = NULL,
    weight.by      = c("prop.var", "eigenvalue", "stdev", "none"),
    k.param        = 20L,
    reduction.name = "wt.pca",
    reduction.key  = "wtPCA_",
    graph.name     = "wt",
    weight.factor  = 1,
    verbose        = TRUE,
    ...
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop(
      "The 'Seurat' package is required for RunWeightedNeighbors(). ",
      "Install it with:\n  install.packages('Seurat')",
      call. = FALSE
    )
  }

  weight.by <- match.arg(weight.by)

  if (!is.numeric(weight.factor) || length(weight.factor) != 1 ||
      weight.factor < 0 || weight.factor > 1) {
    stop("'weight.factor' must be a single number between 0 and 1.",
         call. = FALSE)
  }

  # ── Compute weighted embeddings ────────────────────────────────────────────
  wt <- .compute_weighted_embeddings(object, reduction, dims, weight.by,
                                      weight.factor, verbose)

  # ── Store weighted embeddings as a new DimReduc ───────────────────────────
  source_assay <- tryCatch(
    DefaultAssay(object[[reduction]]),
    error = function(e) DefaultAssay(object)
  )

  wt_emb <- wt$weighted_emb
  rownames(wt_emb) <- rownames(wt$emb)
  colnames(wt_emb) <- paste0(reduction.key, seq_len(ncol(wt_emb)))

  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = wt_emb,
    key        = reduction.key,
    assay      = source_assay,
    misc       = list(
      weight.by        = weight.by,
      weight.factor    = weight.factor,
      weights          = wt$weights,
      source.reduction = reduction,
      dims.used        = wt$dims
    )
  )

  if (verbose) {
    message(sprintf("[wUMAP] Stored weighted embeddings as '%s'", reduction.name))
    message(sprintf("[wUMAP] Running FindNeighbors (k = %d) on '%s'",
                    k.param, reduction.name))
  }

  # ── Build KNN / SNN graphs via Seurat ────────────────────────────────────
  object <- Seurat::FindNeighbors(
    object,
    reduction  = reduction.name,
    dims       = seq_len(ncol(wt_emb)),
    k.param    = k.param,
    graph.name = c(paste0(graph.name, "_nn"), paste0(graph.name, "_snn")),
    verbose    = verbose,
    ...
  )

  object
}
