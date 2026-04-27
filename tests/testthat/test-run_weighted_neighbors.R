# Shared fixture ---------------------------------------------------------------
.make_pbmc <- function(npcs = 5) {
  suppressPackageStartupMessages(library(Seurat))
  data("pbmc_small", package = "SeuratObject")
  suppressWarnings(RunPCA(pbmc_small, npcs = npcs, verbose = FALSE))
}

# RunWeightedNeighbors ---------------------------------------------------------

test_that("RunWeightedNeighbors stores weighted DimReduc and two graphs", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  pbmc <- .make_pbmc()

  result <- suppressWarnings(
    RunWeightedNeighbors(pbmc, dims = 1:5, weight.by = "pct.var",
                         k.param = 5L, graph.name = "wt", verbose = FALSE)
  )

  # Weighted embedding stored
  expect_true("wt.pca" %in% names(result@reductions))
  wt_emb <- Embeddings(result[["wt.pca"]])
  expect_equal(nrow(wt_emb), ncol(pbmc))
  expect_equal(ncol(wt_emb), 5L)

  # Both graphs present
  expect_true("wt_nn"  %in% names(result@graphs))
  expect_true("wt_snn" %in% names(result@graphs))
})

test_that("RunWeightedNeighbors misc slot is populated", {
  skip_if_not_installed("Seurat")

  pbmc <- .make_pbmc()

  result <- suppressWarnings(
    RunWeightedNeighbors(pbmc, dims = 1:5, weight.by = "eigenvalue",
                         k.param = 5L, graph.name = "wt", verbose = FALSE)
  )

  misc <- Misc(result[["wt.pca"]])
  expect_equal(misc$weight.by, "eigenvalue")
  expect_length(misc$weights, 5L)
  expect_equal(misc$source.reduction, "pca")
  expect_equal(misc$dims.used, 1:5)
})

test_that("RunWeightedNeighbors errors on missing reduction", {
  skip_if_not_installed("Seurat")

  pbmc <- .make_pbmc()

  expect_error(
    RunWeightedNeighbors(pbmc, reduction = "nonexistent", verbose = FALSE),
    regexp = "not found"
  )
})

# RunWeightedUMAP — graph path -------------------------------------------------

test_that("RunWeightedUMAP graph path produces correct embedding shape", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  pbmc <- .make_pbmc()

  pbmc <- suppressWarnings(
    RunWeightedNeighbors(pbmc, dims = 1:5, weight.by = "pct.var",
                         k.param = 5L, graph.name = "wt", verbose = FALSE)
  )

  result <- RunWeightedUMAP(
    pbmc,
    graph          = "wt_nn",
    n.neighbors    = 5L,
    n.components   = 2L,
    reduction.name = "wt.umap",
    verbose        = FALSE
  )

  expect_true("wt.umap" %in% names(result@reductions))
  emb <- Embeddings(result[["wt.umap"]])
  expect_equal(nrow(emb), ncol(pbmc))
  expect_equal(ncol(emb), 2L)
  expect_true(all(startsWith(colnames(emb), "wtUMAP_")))
})

test_that("RunWeightedUMAP graph path errors on missing graph", {
  skip_if_not_installed("Seurat")

  pbmc <- .make_pbmc()

  expect_error(
    RunWeightedUMAP(pbmc, graph = "nonexistent_nn", verbose = FALSE),
    regexp = "not found"
  )
})

test_that("graph path and standard path produce same embedding dimensions", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  pbmc <- .make_pbmc()

  pbmc <- suppressWarnings(
    RunWeightedNeighbors(pbmc, dims = 1:5, weight.by = "pct.var",
                         k.param = 5L, graph.name = "wt", verbose = FALSE)
  )

  via_graph <- RunWeightedUMAP(pbmc, graph = "wt_nn", n.neighbors = 5L,
                                n.components = 2L, reduction.name = "umap.graph",
                                verbose = FALSE)

  via_emb   <- RunWeightedUMAP(pbmc, dims = 1:5, weight.by = "pct.var",
                                n.neighbors = 5L, n.components = 2L,
                                reduction.name = "umap.emb", verbose = FALSE)

  expect_equal(dim(Embeddings(via_graph[["umap.graph"]])),
               dim(Embeddings(via_emb[["umap.emb"]])))
})
