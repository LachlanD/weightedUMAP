test_that("RunLocalPCAUMAP returns a Seurat object with correct embedding shape", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- RunLocalPCAUMAP(
    pbmc_small,
    dims         = 1:5,
    k.param      = 10L,
    n.components = 2L,
    verbose      = FALSE
  )

  expect_true("lp.umap" %in% names(result@reductions))

  emb <- Embeddings(result[["lp.umap"]])
  expect_equal(nrow(emb), ncol(pbmc_small))
  expect_equal(ncol(emb), 2L)
  expect_true(all(startsWith(colnames(emb), "lpUMAP_")))
})

test_that("RunLocalPCAUMAP custom reduction.name is stored correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- RunLocalPCAUMAP(
    pbmc_small,
    dims           = 1:5,
    k.param        = 10L,
    reduction.name = "my.lpumap",
    reduction.key  = "myLP_",
    verbose        = FALSE
  )

  expect_true("my.lpumap" %in% names(result@reductions))
  emb <- Embeddings(result[["my.lpumap"]])
  expect_true(all(startsWith(colnames(emb), "myLP_")))
})

test_that("RunLocalPCAUMAP errors on invalid k.param", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  expect_error(
    RunLocalPCAUMAP(pbmc_small, dims = 1:5, k.param = 1L, verbose = FALSE),
    regexp = "k.param"
  )
})

test_that("RunLocalPCAUMAP errors on invalid local.dims", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  expect_error(
    RunLocalPCAUMAP(pbmc_small, dims = 1:5, k.param = 10L,
                    local.dims = 999L, verbose = FALSE),
    regexp = "local.dims"
  )
})

test_that("RunLocalPCAUMAP errors on missing reduction", {
  skip_if_not_installed("Seurat")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")

  expect_error(
    RunLocalPCAUMAP(pbmc_small, reduction = "nonexistent", verbose = FALSE),
    regexp = "not found"
  )
})

test_that("RunLocalPCAUMAP misc slot stores parameters", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- RunLocalPCAUMAP(
    pbmc_small,
    dims       = 1:5,
    k.param    = 10L,
    local.dims = 3L,
    verbose    = FALSE
  )

  misc <- Misc(result[["lp.umap"]])
  expect_equal(misc$k.param, 10L)
  expect_equal(misc$local.dims, 3L)
  expect_equal(misc$source.reduction, "pca")
  expect_equal(misc$dims.used, 1:5)
})
