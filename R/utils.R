# Internal helpers — not exported

# ── .compute_weighted_embeddings ──────────────────────────────────────────────
# Validates inputs, extracts PC embeddings from a Seurat reduction, and applies
# variance-derived weights.  Returns a named list:
#   emb          – raw embeddings (cells × dims)
#   weighted_emb – scaled embeddings (cells × dims)
#   weights      – numeric vector, one weight per dim
#   dims         – integer vector of dim indices used
#   sdev         – standard deviations of the used dims
#
# @keywords internal
.compute_weighted_embeddings <- function(object, reduction, dims, weight.by,
                                         weight.factor, verbose) {
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

  emb_full  <- Embeddings(object[[reduction]])
  sdev_full <- Stdev(object[[reduction]])

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

  weights <- switch(
    weight.by,
    prop.var   = sdev^2 / sum(sdev^2),
    eigenvalue = sdev^2,
    stdev      = sdev,
    none       = rep(1.0, length(sdev))
  )

  # ── Apply weight.factor interpolation ──────────────────────────────────────
  # At weight.factor = 0: all PCs contribute equally (standard UMAP).
  # At weight.factor = 1: full variance-derived weights applied.
  if (weight.factor < 1) {
    weights <- (1 - weight.factor) * mean(weights) + weight.factor * weights
  }

  if (verbose) {
    message(sprintf("[wUMAP] Weight scheme : %s (factor = %.2g)",
                    weight.by, weight.factor))
    message(sprintf("[wUMAP] Dimensions    : %d (dims %d\u2013%d)",
                    length(dims), min(dims), max(dims)))
    max_show <- min(length(weights), 10L)
    wt_str   <- paste(
      sprintf("PC%d=%.3g", dims[seq_len(max_show)], weights[seq_len(max_show)]),
      collapse = ", "
    )
    if (length(weights) > 10L) wt_str <- paste0(wt_str, ", ...")
    message(sprintf("[wUMAP] Weights       : %s", wt_str))
  }

  list(
    emb          = emb,
    weighted_emb = sweep(emb, 2L, weights, `*`),
    weights      = weights,
    dims         = dims,
    sdev         = sdev
  )
}
