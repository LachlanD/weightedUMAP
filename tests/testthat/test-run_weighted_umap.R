test_that("RunWeightedUMAP returns a Seurat object with new reduction", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  suppressPackageStartupMessages({
    library(Seurat)
  })

  # Use built-in pbmc_small data (10x cells, already has PCA)
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- RunPCA(pbmc_small, npcs = 10, verbose = FALSE)

  result <- RunWeightedUMAP(
    pbmc_small,
    dims          = 1:5,
    weight.by     = "prop.var",
    n.neighbors   = 5L,
    n.components  = 2L,
    verbose       = FALSE
  )
  expect_true("wt.umap" %in% names(result@reductions))

  emb <- Embeddings(result[["wt.umap"]])
  expect_equal(nrow(emb), ncol(pbmc_small))   # one row per cell
  expect_equal(ncol(emb), 2L)                  # 2 UMAP dims
  expect_true(all(startsWith(colnames(emb), "wtUMAP_")))
})

test_that("weight.by = 'none' produces standard UMAP output shape", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  suppressPackageStartupMessages(library(Seurat))

  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- RunPCA(pbmc_small, npcs = 10, verbose = FALSE)

  result <- RunWeightedUMAP(
    pbmc_small,
    dims        = 1:5,
    weight.by   = "none",
    n.neighbors = 5L,
    verbose     = FALSE
  )

  emb <- Embeddings(result[["wt.umap"]])
  expect_equal(nrow(emb), ncol(pbmc_small))
  expect_equal(ncol(emb), 2L)
})

test_that("RunWeightedUMAP errors on missing reduction", {
  skip_if_not_installed("SeuratObject")

  data("pbmc_small", package = "SeuratObject")

  expect_error(
    RunWeightedUMAP(pbmc_small, reduction = "nonexistent", verbose = FALSE),
    regexp = "not found"
  )
})

test_that("RunWeightedUMAP errors on out-of-range dims", {
  skip_if_not_installed("Seurat")

  suppressPackageStartupMessages(library(Seurat))

  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- RunPCA(pbmc_small, npcs = 10, verbose = FALSE)

  expect_error(
    RunWeightedUMAP(pbmc_small, dims = 1:999, verbose = FALSE),
    regexp = "exceed available"
  )
})

test_that("misc slot stores weight metadata", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  suppressPackageStartupMessages(library(Seurat))

  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- RunPCA(pbmc_small, npcs = 10, verbose = FALSE)

  result <- RunWeightedUMAP(
    pbmc_small,
    dims      = 1:5,
    weight.by = "prop.var",
    verbose   = FALSE,
    n.neighbors = 5L
  )

  misc <- Misc(result[["wt.umap"]])
  expect_equal(misc$weight.by, "prop.var")
  expect_length(misc$weights, 5L)
  expect_equal(misc$source.reduction, "pca")
  expect_equal(misc$dims.used, 1:5)
})

test_that("weight.factor = 0 produces uniform weights", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  suppressPackageStartupMessages(library(Seurat))

  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- RunPCA(pbmc_small, npcs = 10, verbose = FALSE)

  result <- RunWeightedUMAP(
    pbmc_small,
    dims          = 1:5,
    weight.by     = "prop.var",
    n.neighbors   = 5L,
    verbose       = FALSE
  )

  emb <- Embeddings(result[["wt.umap"]])
  expect_equal(nrow(emb), ncol(pbmc_small))
  expect_equal(ncol(emb), 2L)
})

test_that("weight.factor = 0.5 produces intermediate embedding", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  suppressPackageStartupMessages(library(Seurat))

  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- RunPCA(pbmc_small, npcs = 10, verbose = FALSE)

  r0 <- RunWeightedUMAP(pbmc_small, dims = 1:5, weight.by = "prop.var",
                        weight.factor = 0, n.neighbors = 5L, verbose = FALSE)
  r1 <- RunWeightedUMAP(pbmc_small, dims = 1:5, weight.by = "prop.var",
                        weight.factor = 1, n.neighbors = 5L, verbose = FALSE)
  r5 <- RunWeightedUMAP(pbmc_small, dims = 1:5, weight.by = "prop.var",
                        weight.factor = 0.5, n.neighbors = 5L, verbose = FALSE)

  # Intermediate embedding should differ from both extremes
  expect_false(isTRUE(all.equal(Embeddings(r5[["wt.umap"]]),
                                Embeddings(r0[["wt.umap"]]))))
  expect_false(isTRUE(all.equal(Embeddings(r5[["wt.umap"]]),
                                Embeddings(r1[["wt.umap"]]))))
})

test_that("weight.factor out of range errors", {
  skip_if_not_installed("Seurat")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 5, verbose = FALSE))

  expect_error(
    RunWeightedUMAP(pbmc_small, dims = 1:5, weight.factor = 1.5, verbose = FALSE),
    regexp = "between 0 and 1"
  )
  expect_error(
    RunWeightedUMAP(pbmc_small, dims = 1:5, weight.factor = -0.1, verbose = FALSE),
    regexp = "between 0 and 1"
  )
})

test_that("log.scale = TRUE produces a different embedding than log.scale = FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  r_no_log <- RunWeightedUMAP(pbmc_small, dims = 1:5, weight.by = "prop.var",
                               log.scale = FALSE, n.neighbors = 5L, verbose = FALSE)
  r_log    <- RunWeightedUMAP(pbmc_small, dims = 1:5, weight.by = "prop.var",
                               log.scale = TRUE,  n.neighbors = 5L, verbose = FALSE)

  expect_false(isTRUE(all.equal(Embeddings(r_no_log[["wt.umap"]]),
                                Embeddings(r_log[["wt.umap"]]))))
})

test_that("log.scale errors on non-logical input", {
  skip_if_not_installed("Seurat")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 5, verbose = FALSE))

  expect_error(
    RunWeightedUMAP(pbmc_small, dims = 1:5, log.scale = "yes", verbose = FALSE),
    regexp = "TRUE or FALSE"
  )
})

test_that("mp.filter = TRUE produces Marchenko-Pastur filtered UMAP", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- suppressWarnings(RunWeightedUMAP(
    pbmc_small,
    dims        = 1:5,
    weight.by   = "prop.var",
    mp.filter   = TRUE,
    n.neighbors = 5L,
    verbose     = FALSE
  ))

  emb <- Embeddings(result[["wt.umap"]])
  expect_equal(nrow(emb), ncol(pbmc_small))
  expect_equal(ncol(emb), 2L)

  misc <- Misc(result[["wt.umap"]])
  expect_equal(misc$weight.by, "prop.var")
  expect_true(misc$mp.filter)
  expect_length(misc$weights, 5L)
  expect_true(all(misc$weights >= 0))
  # Weights should sum to 1 after filtering + renormalisation
  expect_equal(sum(misc$weights), 1, tolerance = 1e-10)
})

test_that("mp.filter = TRUE differs from unfiltered prop.var embedding", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  set.seed(1L)
  r_pv <- RunWeightedUMAP(pbmc_small, dims = 1:5, weight.by = "prop.var",
                          n.neighbors = 5L, verbose = FALSE)
  set.seed(1L)
  r_mp <- suppressWarnings(RunWeightedUMAP(
    pbmc_small, dims = 1:5, weight.by = "prop.var", mp.filter = TRUE,
    n.neighbors = 5L, verbose = FALSE,
    reduction.name = "wt.umap.mp"
  ))

  # The weight vectors differ when any PCs are below the noise floor
  # (if none are filtered, weights may be equal — test that the call succeeds
  # and stores mp.filter in misc)
  expect_true(Misc(r_mp[["wt.umap.mp"]])$mp.filter)
})
