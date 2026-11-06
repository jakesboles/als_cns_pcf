dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
library(tidyverse)
library(SpatialExperiment)
library(scater)
library(dittoSeq)
library(pheatmap)

setwd("/projects/b1169/boles/als_cns_pcf")

results_dir <- "results/08_loadings_clustering/"
data_dir <- "data/08_loadings_clustering/"
dir.create(results_dir, showWarnings = F, recursive = T)
dir.create(data_dir, showWarnings = F, recursive = T)

num_pcs <- 12
cluster_on <- "HARMONY"
leiden_res <- c(0.1, 0.3, 0.5, 1)
snn_k <- 15

sce <- readRDS("data/05_object_creation/spe.rds")
sce <- sce[, sce$keep]

rd <- c(readRDS("data/07_embedding/reduced_dims.rds"),
        readRDS("data/07_embedding/reduced_dims_harmony.rds"))
for (nm in names(rd)) {
  if (!is.null(rownames(rd[[nm]]))) stopifnot(identical(rownames(rd[[nm]]), colnames(sce)))
  stopifnot(nrow(rd[[nm]]) == ncol(sce))
  reducedDim(sce, nm) <- rd[[nm]]
}
rm(rd)

# PCA loadings ------------------------------------------------------------

stopifnot(!is.null(attr(reducedDim(sce, "PCA"), "rotation")),
          setequal(rownames(attr(reducedDim(sce, "PCA"), "rotation")), rownames(sce)))
rotation <- attr(reducedDim(sce, "PCA"), "rotation")[, 1:num_pcs]
percent_var <- attr(reducedDim(sce, "PCA"), "percentVar")[1:num_pcs]

write.csv(rotation, paste0(results_dir, "pca_loadings.csv"))

png(paste0(results_dir, "pca_loadings_heatmap.png"), units = "in", res = 300, height = 14, width = 7, bg = "white")
pheatmap(rotation,
         color = colorRampPalette(c("dodgerblue", "white", "firebrick"))(101),
         breaks = seq(-max(abs(rotation)), max(abs(rotation)), length.out = 102),
         cluster_cols = F,
         labels_col = paste0(colnames(rotation), " (", round(percent_var, 1), "%)"),
         fontsize_row = 7)
dev.off()

top_loadings <- as.data.frame(rotation) %>%
  rownames_to_column("marker") %>%
  pivot_longer(-marker, names_to = "pc", values_to = "loading") %>%
  mutate(pc = factor(pc, levels = colnames(rotation))) %>%
  group_by(pc) %>%
  slice_max(abs(loading), n = 10) %>%
  ungroup() %>%
  arrange(pc, loading) %>%
  mutate(y = factor(paste(marker, pc, sep = "___"), levels = unique(paste(marker, pc, sep = "___"))))

p <- top_loadings %>%
  ggplot(aes(x = loading, y = y, fill = loading > 0)) +
  geom_col(show.legend = F) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  facet_wrap(. ~ pc, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "dodgerblue")) +
  labs(y = NULL) +
  theme_linedraw(base_size = 8)
ggsave(paste0(results_dir, "pca_top_loadings.png"), p, width = 12, height = 10, dpi = 300)

# Leiden clustering -------------------------------------------------------

set.seed(1)
graph <- bluster::makeSNNGraph(reducedDim(sce, cluster_on)[, 1:num_pcs],
                               k = snn_k,
                               BNPARAM = BiocNeighbors::AnnoyParam())

leiden <- function(g, res) {
  args <- list(g, objective_function = "modularity", n_iterations = 5)
  if ("resolution" %in% names(formals(igraph::cluster_leiden))) args$resolution <- res else args$resolution_parameter <- res
  set.seed(1)
  factor(igraph::membership(do.call(igraph::cluster_leiden, args)))
}

for (res in leiden_res) {
  col <- paste0("leiden_", res)
  colData(sce)[[col]] <- leiden(graph, res)
  message(col, ": ", nlevels(colData(sce)[[col]]), " clusters")
}
rm(graph)
gc()

cluster_cols <- paste0("leiden_", leiden_res)
arrow::write_parquet(as.data.frame(colData(sce))[, c("cell_id", cluster_cols)],
                     paste0(data_dir, "clusters.parquet"))

# Cluster QC --------------------------------------------------------------

umap_name <- paste0(cluster_on, "_UMAP")

for (col in cluster_cols) {
  cl <- colData(sce)[[col]]

  p <- dittoDimPlot(sce, col, reduction.use = umap_name, do.label = T, labels.size = 3, legend.show = F)
  ggsave(paste0(results_dir, col, "_umap.png"), p, width = 7, height = 6, dpi = 300)

  ct <- table(cluster = cl, cell_type = sce$cell_type)
  write.csv(as.data.frame.matrix(ct), paste0(results_dir, col, "_celltype_counts.csv"))

  png(paste0(results_dir, col, "_celltype_fraction.png"), units = "in", res = 300,
      height = max(4, nlevels(cl) * 0.2), width = 6, bg = "white")
  pheatmap(prop.table(ct, 1), cluster_cols = F, display_numbers = T, number_format = "%.2f", fontsize_number = 6,
           main = "Fraction of each cluster by HALO cell type")
  dev.off()

  means <- sapply(split(seq_len(ncol(sce)), cl), function(i) rowMeans(assay(sce, "asinh")[, i, drop = F]))
  png(paste0(results_dir, col, "_marker_means.png"), units = "in", res = 300,
      height = 12, width = max(5, nlevels(cl) * 0.25 + 3), bg = "white")
  pheatmap(means, scale = "row", fontsize_row = 7, main = "Mean asinh per cluster (row z-score)")
  dev.off()
}
