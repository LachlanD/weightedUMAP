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

# ── local.weight.by tests ──────────────────────────────────────────────────────

test_that("RunLocalPCAUMAP local.weight.by = 'prop.var' stores correct misc", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- RunLocalPCAUMAP(
    pbmc_small,
    dims            = 1:5,
    k.param         = 10L,
    local.dims      = 3L,
    local.weight.by = "prop.var",
    verbose         = FALSE
  )

  emb  <- Embeddings(result[["lp.umap"]])
  misc <- Misc(result[["lp.umap"]])
  expect_equal(nrow(emb), ncol(pbmc_small))
  expect_equal(ncol(emb), 2L)
  expect_equal(misc$local.weight.by, "prop.var")
})

test_that("RunLocalPCAUMAP mp.filter = TRUE produces a valid embedding", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- RunLocalPCAUMAP(
    pbmc_small,
    dims            = 1:5,
    k.param         = 10L,
    local.weight.by = "prop.var",
    mp.filter       = TRUE,
    reduction.name  = "lp.mp",
    verbose         = FALSE
  )

  emb  <- Embeddings(result[["lp.mp"]])
  misc <- Misc(result[["lp.mp"]])
  expect_equal(nrow(emb), ncol(pbmc_small))
  expect_equal(ncol(emb), 2L)
  expect_true(misc$mp.filter)
})

test_that("RunLocalPCAUMAP 'none' and 'prop.var' produce different embeddings", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  set.seed(1L)
  r_none <- RunLocalPCAUMAP(pbmc_small, dims = 1:5, k.param = 10L,
                             local.weight.by = "none",     verbose = FALSE)
  set.seed(1L)
  r_pv   <- RunLocalPCAUMAP(pbmc_small, dims = 1:5, k.param = 10L,
                             local.weight.by = "prop.var", verbose = FALSE,
                             reduction.name = "lp.pv")

  expect_false(isTRUE(all.equal(
    Embeddings(r_none[["lp.umap"]]),
    Embeddings(r_pv[["lp.pv"]])
  )))
})

# ── RunLocalPCANeighbors tests ─────────────────────────────────────────────────

test_that("RunLocalPCANeighbors stores knn and snn graphs", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- RunLocalPCANeighbors(
    pbmc_small,
    dims    = 1:5,
    k.param = 10L,
    prefix  = "lp",
    verbose = FALSE
  )

  expect_true("lp_nn"  %in% names(result@graphs))
  expect_true("lp_snn" %in% names(result@graphs))

  n <- ncol(pbmc_small)
  expect_equal(dim(result[["lp_nn"]]),  c(n, n))
  expect_equal(dim(result[["lp_snn"]]), c(n, n))
})

test_that("RunLocalPCANeighbors knn graph is binary (0/1)", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- RunLocalPCANeighbors(pbmc_small, dims = 1:5, k.param = 10L,
                                  prefix = "lp", verbose = FALSE)

  nn_vals <- result[["lp_nn"]]@x
  expect_true(all(nn_vals %in% c(0L, 1L)))
})

test_that("RunLocalPCANeighbors snn graph has values in [0, 1]", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- RunLocalPCANeighbors(pbmc_small, dims = 1:5, k.param = 10L,
                                  prefix = "lp", verbose = FALSE)

  snn_vals <- result[["lp_snn"]]@x
  expect_true(all(snn_vals >= 0 & snn_vals <= 1))
})

test_that("RunLocalPCANeighbors with local.weight.by = 'prop.var' succeeds", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- RunLocalPCANeighbors(pbmc_small, dims = 1:5, k.param = 10L,
                                  local.weight.by = "prop.var",
                                  prefix = "lp.wt", verbose = FALSE)

  expect_true("lp.wt_nn"  %in% names(result@graphs))
  expect_true("lp.wt_snn" %in% names(result@graphs))
})

test_that("RunLocalPCANeighbors snn graph is usable with FindClusters", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("RANN")

  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  pbmc_small <- suppressWarnings(RunPCA(pbmc_small, npcs = 10, verbose = FALSE))

  result <- RunLocalPCANeighbors(pbmc_small, dims = 1:5, k.param = 10L,
                                  prefix = "lp", verbose = FALSE)

  clustered <- FindClusters(result, graph.name = "lp_snn",
                             resolution = 0.5, verbose = FALSE)
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
  expect_equal(length(clustered$seurat_clusters), ncol(pbmc_small))
})

test_that("RunLocalPCANeighbors errors on invalid reduction", {
  skip_if_not_installed("RANN")

  data("pbmc_small", package = "SeuratObject")

  expect_error(
    RunLocalPCANeighbors(pbmc_small, reduction = "nonexistent"),
    regexp = "not found"
  )
})

