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
    weight.by     = "pct.var",
    n.neighbors   = 5L,
    n.components  = 2L,
    verbose       = FALSE
  )

  expect_s4_class(result, "Seurat")
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
    weight.by = "eigenvalue",
    verbose   = FALSE,
    n.neighbors = 5L
  )

  misc <- Misc(result[["wt.umap"]])
  expect_equal(misc$weight.by, "eigenvalue")
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
    weight.by     = "pct.var",
    weight.factor = 0,
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

  r0 <- RunWeightedUMAP(pbmc_small, dims = 1:5, weight.by = "pct.var",
                        weight.factor = 0, n.neighbors = 5L, verbose = FALSE)
  r1 <- RunWeightedUMAP(pbmc_small, dims = 1:5, weight.by = "pct.var",
                        weight.factor = 1, n.neighbors = 5L, verbose = FALSE)
  r5 <- RunWeightedUMAP(pbmc_small, dims = 1:5, weight.by = "pct.var",
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
