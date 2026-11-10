dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
suppressPackageStartupMessages({
  library(tidyverse)
  library(SpatialExperiment)
  library(scater)
  library(presto)
  library(dittoSeq)
  library(harmony)
  library(pheatmap)
})

setwd("/projects/b1169/boles/als_cns_pcf")

# dataset_z: 05's z (raw, all cells, whole dataset; Hickey 2021); log_section_z: log1p, z per section, capped
norms <- c("dataset_z", "log_section_z")
clust_norm <- Sys.getenv("PCF_NORM", "dataset_z")
if (!clust_norm %in% norms) stop("PCF_NORM must be one of ", paste(norms, collapse = ", "))
message("Clustering input: ", clust_norm)

# panel: pca_features only; all: every channel kept in 05
feature_set <- Sys.getenv("PCF_FEATURES", "panel")
if (!feature_set %in% c("panel", "all")) stop("PCF_FEATURES must be panel or all")
message("PCA features: ", feature_set)

# knn: UMAP-style fuzzy kNN graph (scanpy sc.pp.neighbors, as in the AD paper); snn: bluster shared-neighbour graph
graph_type <- Sys.getenv("PCF_GRAPH", "snn")
if (!graph_type %in% c("knn", "snn")) stop("PCF_GRAPH must be knn or snn")
graph_k <- as.integer(Sys.getenv("PCF_K", if (graph_type == "knn") "150" else "15"))
run_name <- paste0(clust_norm, if (graph_type == "snn" && graph_k == 15) "" else paste0("_", graph_type, graph_k),
                   if (feature_set == "all") "_all_features" else "")
message("Graph: ", graph_type, ", k = ", graph_k, "; outputs: ", run_name)

results_dir <- paste0("results/07b_all_cells/", run_name, "/")
data_dir <- paste0("data/07b_all_cells/", run_name, "/all/")
dir.create(results_dir, showWarnings = F, recursive = T)
dir.create(data_dir, showWarnings = F, recursive = T)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

pca_features <- c("TMEM119", "Olig-2", "CD3e", "CD4",
                  "NeuN", "GFAP", "Vimentin", "CD11c",
                  "Iba-1", "CD163", "CD8")
num_pcs <- as.integer(Sys.getenv("PCF_NPCS", if (feature_set == "all") "12" else as.character(length(pca_features) - 1)))
leiden_res <- c(0.1, 0.3, 0.5, 1)
z_cap <- 10

sce <- readRDS("data/05_object_creation/spe.rds")
message("Removing ", sum(sce$flag_size %in% TRUE), " cells that failed the size filter")
sce <- sce[, !(sce$flag_size %in% TRUE)]
if (clust_norm == "dataset_z" && !"z" %in% assayNames(sce)) stop("no z assay in the SPE; rerun 05 to the end")
sce$cell_type <- factor(sce$cell_type)
sce$code <- factor(sce$code)

missing <- setdiff(pca_features, rownames(sce))
if (length(missing)) stop("pca_features not in rownames(sce): ", paste(missing, collapse = ", "))
clust_features <- if (feature_set == "all") rownames(sce) else pca_features
message(length(clust_features), " PCA features, ", num_pcs, " PCs used")

# clustering input, PCA features only
if (clust_norm == "dataset_z") {
  clust <- as.matrix(assay(sce, "z")[clust_features, , drop = F])
} else {
  clust <- log1p(pmax(as.matrix(assay(sce, "counts")[clust_features, , drop = F]), 0))
  for (cd in levels(sce$code)) {
    i <- which(sce$code == cd)
    mu <- rowMeans(clust[, i, drop = F])
    sdv <- matrixStats::rowSds(clust[, i, drop = F])
    sdv[is.na(sdv) | sdv == 0] <- 1
    clust[, i] <- pmax(pmin((clust[, i, drop = F] - mu) / sdv, z_cap), -z_cap)
  }
}
assays(sce) <- assays(sce)["asinh"]
gc()
sce$flag_status <- factor(case_when(is.na(sce$keep) ~ "NA",
                                    sce$keep ~ "kept",
                                    TRUE ~ sce$flag_reason))
sce
print(table(sce$cell_type, sce$flag_status))

save_umap <- function(p, file, w){
  ggsave(filename = paste0(results_dir, file),
         plot = p,
         units = "in", dpi = 600,
         height = 5, width = w)
}

set.seed(1)
explain_idx <- sample(ncol(sce), min(ncol(sce), 2e5))

explain_vars <- c("code", "batch", "cell_type", "group", "sex", "age", "area_um2",
                  "keep", "flag_size", "flag_region", "flag_lineage")
explain_vars <- explain_vars[map_lgl(explain_vars, ~ length(unique(na.omit(colData(sce)[[.x]][explain_idx]))) > 1)]
for (v in intersect(explain_vars, c("keep", "flag_size", "flag_region", "flag_lineage")))
  colData(sce)[[v]] <- factor(colData(sce)[[v]])

plot_explanatory <- function(dimred, file) {
  ex <- getExplanatoryPCs(sce[, explain_idx], variables = explain_vars, dimred = dimred)
  png(filename = paste0(results_dir, file), units = "in", res = 600, height = 6, width = 7, bg = "white")
  print(plotExplanatoryPCs(ex / 100))
  dev.off()
}

umap_plots <- function(dimred, prefix) {
  files <- c(cell_type = "celltype", batch = "batch", code = "sample", flag_status = "flags")
  for (v in names(files)) {
    p <- dittoDimPlot(sce, v, reduction.use = dimred, order = "randomize")
    save_umap(p, paste0(prefix, "_umap_", files[[v]], ".png"), if (v == "code") 9 else 6)
  }
}

# Dimensional reduction ---------------------------------------------------

reducedDim(sce, "PCA") <- scater::calculatePCA(clust,
                                               scale = F,
                                               ncomponents = min(length(clust_features), 30))
rm(clust)
gc()
stopifnot(num_pcs <= ncol(reducedDim(sce, "PCA")))

# PCA loadings ------------------------------------------------------------

rotation <- attr(reducedDim(sce, "PCA"), "rotation")
percent_var <- attr(reducedDim(sce, "PCA"), "percentVar")

write.csv(rotation, paste0(results_dir, "pca_loadings.csv"))

png(filename = paste0(results_dir, "pca_elbow.png"),
    units = "in", res = 600,
    height = 4, width = 5,
    bg = "white")
plot(percent_var, log="y", xlab="PC", ylab="Variance explained (%)")
dev.off()

png(paste0(results_dir, "pca_loadings_heatmap.png"), units = "in", res = 300, height = max(14, nrow(rotation) * 0.2 + 2),
    width = max(7, ncol(rotation) * 0.35 + 3), bg = "white")
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
  filter(as.integer(pc) <= 20) %>%
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
ggsave(paste0(results_dir, "pca_top_loadings.png"), p, width = 12, height = max(10, 2.5 * ceiling(nlevels(droplevels(top_loadings$pc)) / 4)), dpi = 300)

plot_explanatory("PCA", "explanatory_pcs.png")

set.seed(1)
reducedDim(sce, "UMAP") <- uwot::umap(reducedDim(sce, "PCA")[, 1:num_pcs],
                                      n_neighbors = 15,
                                      n_threads = n_cores,
                                      n_sgd_threads = n_cores,
                                      verbose = T)

saveRDS(list(PCA = reducedDim(sce, "PCA"),
             UMAP = reducedDim(sce, "UMAP")),
        paste0(data_dir, "reduced_dims.rds"))

umap_plots("UMAP", "raw")

# Harmony -----------------------------------------------------------------

sce <- RunHarmony(sce,
                  group.by.vars = "code")

plot_explanatory("HARMONY", "explanatory_harmony.png")

set.seed(1)
reducedDim(sce, "HARMONY_UMAP") <- uwot::umap(reducedDim(sce, "HARMONY")[, 1:num_pcs],
                                              n_neighbors = 15,
                                              n_threads = n_cores,
                                              n_sgd_threads = n_cores)

saveRDS(list(HARMONY = reducedDim(sce, "HARMONY"),
             HARMONY_UMAP = reducedDim(sce, "HARMONY_UMAP")),
        paste0(data_dir, "reduced_dims_harmony.rds"))

umap_plots("HARMONY_UMAP", "harmony")

marker_plots <- map(pca_features, ~ dittoDimPlot(sce[, explain_idx], .x, assay = "asinh", reduction.use = "HARMONY_UMAP",
                                                 order = "increasing", size = 0.3, legend.show = F) +
                      labs(title = .x, x = NULL, y = NULL))
p <- gridExtra::arrangeGrob(grobs = marker_plots, ncol = 4)
ggsave(paste0(results_dir, "harmony_umap_markers.png"), p, width = 16, height = 4 * ceiling(length(pca_features) / 4), dpi = 300)

# Leiden clustering -------------------------------------------------------

set.seed(1)
if (graph_type == "knn") {
  if (!exists("similarity_graph", envir = asNamespace("uwot"))) stop("uwot too old for similarity_graph(); update uwot")
  fg <- uwot::similarity_graph(reducedDim(sce, "HARMONY")[, 1:num_pcs], n_neighbors = graph_k,
                               nn_method = "annoy", n_threads = n_cores, verbose = T)
  graph <- igraph::graph_from_adjacency_matrix(fg, mode = "undirected", weighted = T, diag = F)
  rm(fg)
} else {
  snn_args <- list(reducedDim(sce, "HARMONY")[, 1:num_pcs], k = graph_k, BNPARAM = BiocNeighbors::AnnoyParam())
  if ("num.threads" %in% names(formals(bluster::makeSNNGraph))) snn_args$num.threads <- n_cores else
    snn_args$BPPARAM <- BiocParallel::MulticoreParam(n_cores)
  graph <- do.call(bluster::makeSNNGraph, snn_args)
}
message("Graph: ", igraph::vcount(graph), " cells, ", igraph::ecount(graph), " edges")
gc()

leiden <- function(g, res) {
  args <- list(g, objective_function = "modularity", n_iterations = 5)
  if ("resolution" %in% names(formals(igraph::cluster_leiden))) args$resolution <- res else args$resolution_parameter <- res
  set.seed(1)
  factor(igraph::membership(do.call(igraph::cluster_leiden, args)))
}

cluster_cols <- paste0("leiden_", leiden_res)
for (i in seq_along(leiden_res)) {
  colData(sce)[[cluster_cols[i]]] <- leiden(graph, leiden_res[i])
  message(cluster_cols[i], ": ", nlevels(colData(sce)[[cluster_cols[i]]]), " clusters")
}
rm(graph)
gc()

arrow::write_parquet(as.data.frame(colData(sce))[, c("cell_id", "cell_type", cluster_cols)],
                     paste0(data_dir, "clusters.parquet"))

# Cluster QC --------------------------------------------------------------

frac_heatmap <- function(tab, file, main) {
  png(paste0(results_dir, file), units = "in", res = 300,
      height = max(4, nrow(tab) * 0.2 + 1.5), width = max(6, ncol(tab) * 0.3 + 2), bg = "white")
  pheatmap(prop.table(tab, 1), cluster_rows = F, cluster_cols = F, display_numbers = T,
           number_format = "%.2f", fontsize_number = 5, main = main)
  dev.off()
}

for (col in cluster_cols) {
  cl <- colData(sce)[[col]]

  p <- dittoDimPlot(sce, col, reduction.use = "HARMONY_UMAP", do.label = T, labels.size = 3,
                    legend.show = F, order = "randomize")
  ggsave(paste0(results_dir, col, "_umap.png"), p, width = 7, height = 6, dpi = 300)

  for (v in c("code", "flag_status", "cell_type")) {
    ct <- table(cluster = cl, colData(sce)[[v]])
    write.csv(as.data.frame.matrix(ct), paste0(results_dir, col, "_", v, "_counts.csv"))
    frac_heatmap(ct, paste0(col, "_", v, "_fraction.png"), paste("Fraction of each cluster by", v))
  }
  frac_heatmap(t(table(cluster = cl, sce$cell_type)), paste0(col, "_cell_type_by_cluster_fraction.png"),
               "Fraction of each cell type in each cluster")

  if (nlevels(cl) < 2) next
  wx <- presto::wilcoxauc(assay(sce, "asinh"), as.character(cl))
  wx$log2FC <- wx$logFC / log(2)
  write.csv(wx, paste0(results_dir, col, "_markers.csv"), row.names = F)
  auc_m <- with(wx, tapply(auc, list(feature, group), identity))[, levels(cl), drop = F]
  lfc_m <- with(wx, tapply(log2FC, list(feature, group), identity))[rownames(auc_m), levels(cl), drop = F]
  ord <- hclust(dist(auc_m))$order
  lim <- max(quantile(abs(lfc_m), 0.99, na.rm = T), 0.1)
  pal <- colorRampPalette(c("dodgerblue", "white", "firebrick"))(101)
  h_auc <- pheatmap(auc_m[ord, ], cluster_rows = F, cluster_cols = F, color = pal,
                    breaks = seq(0, 1, length.out = 102), fontsize_row = 7,
                    main = "AUROC, cluster vs all other cells", silent = T)
  h_lfc <- pheatmap(pmax(pmin(lfc_m[ord, ], lim), -lim), cluster_rows = F, cluster_cols = F, color = pal,
                    breaks = seq(-lim, lim, length.out = 102), fontsize_row = 7,
                    main = sprintf("log2FC vs rest (mean asinh diff / ln 2; capped at %.2f)", lim), silent = T)
  png(paste0(results_dir, col, "_markers_auc_log2fc.png"), units = "in", res = 300,
      height = 12, width = 2 * max(5, nlevels(cl) * 0.3 + 3), bg = "white")
  gridExtra::grid.arrange(h_auc$gtable, h_lfc$gtable, ncol = 2)
  dev.off()
}
