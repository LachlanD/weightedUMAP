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
                                         weight.factor, log.scale, verbose,
                                         mp.filter = FALSE) {
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
    stdev      = sdev / sum(sdev),
    none       = rep(1.0, length(sdev))
  )

  # ── Optional MP filtering ────────────────────────────────────────────────
  # When mp.filter = TRUE, zero out PCs whose variance is at or below the
  # Marchenko-Pastur bulk-noise upper bound (lambda_max), then renormalise.
  if (mp.filter) {
    n <- nrow(emb)
    p <- tryCatch(nrow(object[[reduction]]@feature.loadings),
                  error = function(e) NULL)
    if (is.null(p) || p == 0L) {
      warning(
        "mp.filter = TRUE but could not determine the number of PCA features. ",
        "Filter has no effect.",
        call. = FALSE
      )
    } else {
      gamma   <- p / n
      lam_max <- (1.0 + sqrt(gamma))^2
      signal  <- sdev^2 > lam_max
      if (!any(signal)) {
        warning(
          "mp.filter = TRUE but no PCs exceed the noise threshold ",
          sprintf("(lambda_max = %.3g, n = %d, p = %d). ", lam_max, n, p),
          "Filter has no effect.",
          call. = FALSE
        )
      } else {
        weights[!signal] <- 0
        if (!identical(weight.by, "none")) weights <- weights / sum(weights)
        if (verbose) {
          n_sig <- sum(signal)
          message(sprintf(
            "[wUMAP] MP filter     : lambda_max = %.3g (%d/%d PCs above noise floor, %d zeroed)",
            lam_max, n_sig, length(sdev), sum(!signal)
          ))
        }
      }
    }
  }

  # ── Apply log scaling ────────────────────────────────────────────────────────────────────
  if (log.scale && !identical(weight.by, "none")) {
    weights <- log1p(weights)
    weights <- weights / sum(weights)
  }

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

# ── .compute_local_pca_distances ─────────────────────────────────────────────
# Finds k nearest neighbours in PCA space (RANN), then for each cell fits a
# compact SVD to its neighbourhood and measures the displacement to each
# neighbour in the local PCA basis.  Optionally weights each local PC
# direction by its local variance contribution (prop.var / stdev / mp). This
# is analogous to how .compute_weighted_embeddings weights global PCs.
#
# @param emb          Numeric matrix: n cells × d dimensions.
# @param k.param      Number of nearest neighbours (integer, >= 2).
# @param local.dims   Number of local PC directions to retain (integer).
# @param local.weight.by  One of "none", "prop.var", "stdev", "mp".
# @param verbose      Logical. Print progress messages.
#
# @return Named list:
#   nn_idx     – integer matrix n × k.param (neighbour indices, self excluded)
#   local_dist – numeric matrix n × k.param (local PCA distances)
#
# @importFrom RANN nn2
# @keywords internal
.compute_local_pca_distances <- function(emb, k.param, local.dims,
                                          local.weight.by = "none",
                                          mp.filter = FALSE,
                                          verbose = FALSE) {
  n       <- nrow(emb)
  d       <- ncol(emb)
  n_local <- k.param + 1L        # neighbourhood size including focal cell

  # ── Step 1: global k-NN ──────────────────────────────────────────────────
  if (verbose) message("[wUMAP] Step 1: finding global k-NN ...")
  nn_res <- RANN::nn2(emb, k = n_local)
  nn_idx <- nn_res$nn.idx[, -1L, drop = FALSE]   # n × k, exclude self

  # ── Step 2: local PCA distances ──────────────────────────────────────────
  if (verbose) message("[wUMAP] Step 2: computing local PCA distances ...")
  local_dist <- matrix(0.0, nrow = n, ncol = k.param)

  for (i in seq_len(n)) {
    nbr_idx <- nn_idx[i, ]

    # Neighbourhood matrix: focal cell + its k neighbours, centred
    X_nbr <- emb[c(i, nbr_idx), , drop = FALSE]
    X_ctr <- X_nbr - matrix(colMeans(X_nbr), nrow = nrow(X_nbr),
                             ncol = d, byrow = TRUE)

    # Compact SVD — only right singular vectors (local basis) needed
    sv <- tryCatch(
      svd(X_ctr, nu = 0L, nv = local.dims),
      error = function(e) NULL
    )

    if (is.null(sv)) {
      # Fallback: raw Euclidean distance from focal cell to neighbours
      local_dist[i, ] <- sqrt(rowSums(
        (emb[nbr_idx, , drop = FALSE] -
           matrix(emb[i, ], nrow = k.param, ncol = d, byrow = TRUE))^2
      ))
      next
    }

    V     <- sv$v                                 # d × local.dims  local basis
    xi    <- emb[i, , drop = FALSE]              # 1 × d focal cell
    delta <- emb[nbr_idx, , drop = FALSE] -       # k × d displacements
      matrix(xi, nrow = k.param, ncol = d, byrow = TRUE)
    proj  <- delta %*% V                          # k × local.dims projections

    # ── Apply local PC weighting (analogous to global weight.by) ──────────
    if (local.weight.by != "none" || mp.filter) {
      # Variance explained by each local PC direction
      # svd of X_ctr (nrow = k+1): eigenvalue = sv$d^2 / k
      local_var <- sv$d[seq_len(local.dims)]^2 / max(1L, n_local - 1L)

      if (local.weight.by != "none") {
        w_local <- switch(
          local.weight.by,
          prop.var = {
            tot <- sum(local_var)
            if (tot > 0) local_var / tot else rep(1.0 / local.dims, local.dims)
          },
          stdev = {
            s   <- sqrt(pmax(0.0, local_var))
            tot <- sum(s)
            if (tot > 0) s / tot else rep(1.0 / local.dims, local.dims)
          }
        )
      } else {
        w_local <- rep(1.0 / local.dims, local.dims)
      }

      # ── Optional MP filtering (applied after base weights) ────────────────
      if (mp.filter) {
        gamma   <- d / n_local
        lam_max <- (1.0 + sqrt(gamma))^2
        signal  <- local_var > lam_max
        if (any(signal)) {
          w_local[!signal] <- 0
          if (sum(w_local) > 0)
            w_local <- w_local / sum(w_local)
        }
        # If no signal PCs, fall through with unmodified weights
      }

      # Scale each local PC direction by its weight (mirrors global sweep())
      proj <- proj * matrix(w_local, nrow = k.param,
                            ncol = local.dims, byrow = TRUE)
    }

    local_dist[i, ] <- sqrt(rowSums(proj^2))
  }

  list(nn_idx = nn_idx, local_dist = local_dist)
}
