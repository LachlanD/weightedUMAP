# Benchmarking wUMAP Against Standard UMAP

## Strategy

To evaluate whether variance-weighting improves clustering quality, we
need an **independent ground truth** — one not derived from the RNA
clusters themselves.

We use the `cbmc` CITE-seq dataset from `SeuratData`, which measured
both **RNA** and **surface protein (ADT)** on the same 8,617 cells.
ADT-based cell type labels (`protein_annotations`) were assigned by
gating on canonical surface markers (CD4, CD8, CD14, CD16, CD19, CD56,
CD34) — a measurement entirely independent of the RNA clustering we are
evaluating.

Two metrics are computed for each method:

| Metric | What it measures |
|----|----|
| **ARI** (Adjusted Rand Index) | How well RNA clusters match protein-derived cell types. Range 0–1; higher = better. |
| **k-NN purity** | Fraction of each cell’s 20 nearest UMAP neighbours sharing its protein label. Range 0–1; higher = better. |

## Setup

``` r

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(wUMAP)
library(mclust)   # adjustedRandIndex
library(RANN)     # fast k-NN
```

``` r

data(cbmc)
cbmc <- UpdateSeuratObject(cbmc)

# Keep only cells with protein labels; remove doublets
cbmc <- cbmc[, !cbmc$protein_annotations %in% "T/Mono doublets" &
               !is.na(cbmc$protein_annotations)]

cbmc <- NormalizeData(cbmc, verbose = FALSE)
cbmc <- FindVariableFeatures(cbmc, nfeatures = 2000, verbose = FALSE)
cbmc <- ScaleData(cbmc, verbose = FALSE)
cbmc <- RunPCA(cbmc, npcs = 30, verbose = FALSE)
```

``` r

knn_purity <- function(emb, labels, k = 20) {
  nn  <- RANN::nn2(emb, k = k + 1)$nn.idx[, -1]
  lab <- as.integer(factor(labels))
  mean(sapply(seq_len(nrow(nn)), function(i) mean(lab[nn[i, ]] == lab[i])))
}
```

## Standard UMAP (baseline)

``` r

set.seed(42)
cbmc <- FindNeighbors(cbmc, dims = 1:30, verbose = FALSE)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:30, reduction.name = "umap.std", verbose = FALSE)

ari_std <- adjustedRandIndex(cbmc$seurat_clusters, cbmc$protein_annotations)
pur_std <- knn_purity(Embeddings(cbmc, "umap.std"), cbmc$protein_annotations)
```

## Weighted UMAP variants

``` r

set.seed(42)

# stdev weighting
cbmc <- RunWeightedNeighbors(cbmc, dims = 1:30, weight.by = "stdev",
                             prefix = "wt.sd", verbose = FALSE)
cbmc <- FindClusters(cbmc, graph.name = "wt.sd_snn", resolution = 0.8,
                     verbose = FALSE)
cbmc <- RunWeightedUMAP(cbmc, dims = 1:30, weight.by = "stdev",
                        reduction.name = "umap.sd", verbose = FALSE)
ari_sd <- adjustedRandIndex(cbmc$seurat_clusters, cbmc$protein_annotations)
pur_sd <- knn_purity(Embeddings(cbmc, "umap.sd"), cbmc$protein_annotations)

# stdev + log.scale
cbmc <- RunWeightedNeighbors(cbmc, dims = 1:30, weight.by = "stdev",
                             log.scale = TRUE, prefix = "wt.sd.log",
                             verbose = FALSE)
cbmc <- FindClusters(cbmc, graph.name = "wt.sd.log_snn", resolution = 0.8,
                     verbose = FALSE)
cbmc <- RunWeightedUMAP(cbmc, dims = 1:30, weight.by = "stdev",
                        log.scale = TRUE, reduction.name = "umap.sd.log",
                        verbose = FALSE)
ari_sd_log <- adjustedRandIndex(cbmc$seurat_clusters, cbmc$protein_annotations)
pur_sd_log <- knn_purity(Embeddings(cbmc, "umap.sd.log"), cbmc$protein_annotations)

# prop.var, weight.factor = 0.5
cbmc <- RunWeightedNeighbors(cbmc, dims = 1:30, weight.by = "prop.var",
                             weight.factor = 0.5, prefix = "wt.pv05",
                             verbose = FALSE)
cbmc <- FindClusters(cbmc, graph.name = "wt.pv05_snn", resolution = 0.8,
                     verbose = FALSE)
cbmc <- RunWeightedUMAP(cbmc, dims = 1:30, weight.by = "prop.var",
                        weight.factor = 0.5, reduction.name = "umap.pv05",
                        verbose = FALSE)
ari_pv05 <- adjustedRandIndex(cbmc$seurat_clusters, cbmc$protein_annotations)
pur_pv05 <- knn_purity(Embeddings(cbmc, "umap.pv05"), cbmc$protein_annotations)

# prop.var, weight.factor = 1
cbmc <- RunWeightedNeighbors(cbmc, dims = 1:30, weight.by = "prop.var",
                             prefix = "wt.pv", verbose = FALSE)
cbmc <- FindClusters(cbmc, graph.name = "wt.pv_snn", resolution = 0.8,
                     verbose = FALSE)
cbmc <- RunWeightedUMAP(cbmc, dims = 1:30, weight.by = "prop.var",
                        reduction.name = "umap.pv", verbose = FALSE)
ari_pv <- adjustedRandIndex(cbmc$seurat_clusters, cbmc$protein_annotations)
pur_pv <- knn_purity(Embeddings(cbmc, "umap.pv"), cbmc$protein_annotations)

# prop.var + log.scale
cbmc <- RunWeightedNeighbors(cbmc, dims = 1:30, weight.by = "prop.var",
                             log.scale = TRUE, prefix = "wt.pv.log",
                             verbose = FALSE)
cbmc <- FindClusters(cbmc, graph.name = "wt.pv.log_snn", resolution = 0.8,
                     verbose = FALSE)
cbmc <- RunWeightedUMAP(cbmc, dims = 1:30, weight.by = "prop.var",
                        log.scale = TRUE, reduction.name = "umap.pv.log",
                        verbose = FALSE)
ari_pv_log <- adjustedRandIndex(cbmc$seurat_clusters, cbmc$protein_annotations)
pur_pv_log <- knn_purity(Embeddings(cbmc, "umap.pv.log"), cbmc$protein_annotations)

# Marchenko-Pastur weighting
cbmc <- RunWeightedNeighbors(cbmc, dims = 1:30, weight.by = "mp",
                             prefix = "wt.mp", verbose = FALSE)
cbmc <- FindClusters(cbmc, graph.name = "wt.mp_snn", resolution = 0.8,
                     verbose = FALSE)
cbmc <- RunWeightedUMAP(cbmc, dims = 1:30, weight.by = "mp",
                        reduction.name = "umap.mp", verbose = FALSE)
ari_mp <- adjustedRandIndex(cbmc$seurat_clusters, cbmc$protein_annotations)
pur_mp <- knn_purity(Embeddings(cbmc, "umap.mp"), cbmc$protein_annotations)
```

## Local PCA UMAP

[`RunLocalPCAUMAP()`](https://lachland.github.io/weightedUMAP/reference/RunLocalPCAUMAP.md)
is an embedding-only method — it builds a per-neighbourhood PCA basis
and uses it to compute anisotropic local distances before running UMAP.
There is no companion graph-building step, so only k-NN purity (a direct
measure of UMAP embedding quality) is reported; ARI is left as `NA`.

``` r

set.seed(42)
cbmc <- RunLocalPCAUMAP(cbmc, dims = 1:30, k.param = 30,
                        reduction.name = "umap.lp", verbose = FALSE)
pur_lp <- knn_purity(Embeddings(cbmc, "umap.lp"), cbmc$protein_annotations)
```

## Results

``` r

results <- data.frame(
  Method = c(
    "standard (none)",
    "stdev",
    "stdev + log.scale",
    "prop.var  wf=0.5",
    "prop.var  wf=1",
    "prop.var + log.scale",
    "mp (Marchenko-Pastur)",
    "local PCA UMAP"
  ),
  ARI       = c(ari_std, ari_sd, ari_sd_log, ari_pv05, ari_pv, ari_pv_log,
                ari_mp, NA_real_),
  kNN_purity = c(pur_std, pur_sd, pur_sd_log, pur_pv05, pur_pv, pur_pv_log,
                 pur_mp, pur_lp)
)
knitr::kable(results, digits = 3,
  caption = "Clustering quality vs protein labels on cbmc CITE-seq (n ≈ 7,699 cells). ARI = NA for local PCA UMAP (embedding-only method).")
```

| Method                |   ARI | kNN_purity |
|:----------------------|------:|-----------:|
| standard (none)       | 0.548 |      0.927 |
| stdev                 | 0.504 |      0.927 |
| stdev + log.scale     | 0.503 |      0.929 |
| prop.var wf=0.5       | 0.496 |      0.924 |
| prop.var wf=1         | 0.318 |      0.878 |
| prop.var + log.scale  | 0.398 |      0.880 |
| mp (Marchenko-Pastur) | 0.322 |      0.875 |
| local PCA UMAP        |    NA |      0.935 |

Clustering quality vs protein labels on cbmc CITE-seq (n ≈ 7,699 cells).
ARI = NA for local PCA UMAP (embedding-only method). {.table}

``` r

library(tidyr)

res_long <- tidyr::pivot_longer(results, cols = c("ARI", "kNN_purity"),
                                names_to = "Metric", values_to = "Score")
res_long$Method <- factor(res_long$Method, levels = results$Method)

baseline <- data.frame(
  Metric = c("ARI", "kNN_purity"),
  Score  = c(ari_std, pur_std)
)

ggplot(res_long[!is.na(res_long$Score), ],
       aes(x = Method, y = Score, fill = Method)) +
  geom_col(show.legend = FALSE) +
  geom_hline(data = baseline, aes(yintercept = Score),
             linetype = "dashed", colour = "grey40") +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  labs(x = NULL, y = NULL,
       title = "wUMAP benchmark on cbmc CITE-seq",
       subtitle = "Dashed line = standard UMAP baseline") +
  theme_bw(base_size = 12)
```

![ARI and k-NN purity for each method. Dashed line = standard UMAP
baseline. ARI is not shown for local PCA UMAP (embedding-only
method).](benchmark_files/figure-html/results-plot-1.png)

ARI and k-NN purity for each method. Dashed line = standard UMAP
baseline. ARI is not shown for local PCA UMAP (embedding-only method).

## UMAP layouts coloured by protein label

``` r

Idents(cbmc) <- "protein_annotations"

plots <- list(
  DimPlot(cbmc, reduction = "umap.std",   label = TRUE, repel = TRUE) +
    ggtitle(sprintf("standard  ARI=%.3f", ari_std))    + NoLegend(),
  DimPlot(cbmc, reduction = "umap.sd",    label = TRUE, repel = TRUE) +
    ggtitle(sprintf("stdev  ARI=%.3f", ari_sd))        + NoLegend(),
  DimPlot(cbmc, reduction = "umap.sd.log", label = TRUE, repel = TRUE) +
    ggtitle(sprintf("stdev+log  ARI=%.3f", ari_sd_log))+ NoLegend(),
  DimPlot(cbmc, reduction = "umap.pv05",  label = TRUE, repel = TRUE) +
    ggtitle(sprintf("prop.var wf=0.5  ARI=%.3f", ari_pv05)) + NoLegend(),
  DimPlot(cbmc, reduction = "umap.pv",    label = TRUE, repel = TRUE) +
    ggtitle(sprintf("prop.var wf=1  ARI=%.3f", ari_pv)) + NoLegend(),
  DimPlot(cbmc, reduction = "umap.pv.log", label = TRUE, repel = TRUE) +
    ggtitle(sprintf("prop.var+log  ARI=%.3f", ari_pv_log)) + NoLegend(),
  DimPlot(cbmc, reduction = "umap.mp",    label = TRUE, repel = TRUE) +
    ggtitle(sprintf("mp  ARI=%.3f", ari_mp)) + NoLegend(),
  DimPlot(cbmc, reduction = "umap.lp",    label = TRUE, repel = TRUE) +
    ggtitle(sprintf("local PCA UMAP  purity=%.3f", pur_lp)) + NoLegend()
)

wrap_plots(plots, ncol = 2)
```

![All methods coloured by ADT-derived protein
label.](benchmark_files/figure-html/umap-grid-1.png)

All methods coloured by ADT-derived protein label.

## Interpretation

On this PBMC-like dataset, aggressive `prop.var` weighting (`wf=1`)
**reduces** clustering quality (ARI 0.548 → 0.318). PC 1 in scRNA-seq
captures the monocyte–lymphocyte axis so strongly that further
emphasising it collapses finer lymphocyte subtypes together.

Key takeaways:

- **`stdev`** is the safest default — only a marginal ARI change, same
  k-NN purity as standard.
- **`prop.var` with `weight.factor = 0.5`** (partial blend) avoids the
  worst of over-emphasising PC 1.
- **`log.scale = TRUE`** partially mitigates PC 1 dominance and recovers
  some ARI compared to full `prop.var`.
- **`mp` (Marchenko–Pastur)** automatically zeros out noise PCs and
  retains only those with variance above the random-matrix noise
  ceiling. On datasets where most PCs are genuine signal it behaves like
  `prop.var`; on noisier datasets it can sharpen clusters by discarding
  uninformative PCs entirely.
- **Local PCA UMAP** is not evaluated by ARI (it does not build a
  clustering graph), but its k-NN purity reflects how well the embedding
  locally preserves cell-type neighbourhoods. It can expose
  trajectory-like geometry that global metrics miss.
- k-NN purity differences are generally small across all methods,
  indicating UMAP topology is broadly preserved — the main effect is on
  Louvain cluster resolution.

**Recommendation**: start with `weight.by = "stdev"` or
`weight.by = "prop.var", weight.factor = 0.5`, and inspect cluster
marker genes to judge whether the weighting is helping or obscuring
biology in your dataset. Use `weight.by = "mp"` when you want a
principled, data-driven noise threshold. Use
[`RunLocalPCAUMAP()`](https://lachland.github.io/weightedUMAP/reference/RunLocalPCAUMAP.md)
when you suspect the dataset has strong local anisotropic structure
(trajectories, gradients) that a global metric would miss.
