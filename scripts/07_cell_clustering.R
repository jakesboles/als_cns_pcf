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

# one array task per entry, in this order
cell_groups <- list(
  astrocyte = "astrocyte",
  cortical_neuron = "cortical_neuron",
  myeloid = "myeloid",
  oligodendrocyte = "oligodendrocyte",
  tcell = c("cd4_tcell", "cd8_tcell")
)

cell_group <- Sys.getenv("PCF_GROUP")
if (cell_group == "") {
  task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
  if (is.na(task) || task < 1 || task > length(cell_groups))
    stop("Set PCF_GROUP or run as array task 1-", length(cell_groups), ": ", paste(names(cell_groups), collapse = ", "))
  cell_group <- names(cell_groups)[task]
}
if (!cell_group %in% names(cell_groups)) stop("Unknown PCF_GROUP ", cell_group, "; one of ", paste(names(cell_groups), collapse = ", "))
message("Cell group: ", cell_group, " (", paste(cell_groups[[cell_group]], collapse = ", "), ")")

results_dir <- paste0("results/07_cell_clustering/", cell_group, "/")
data_dir <- paste0("data/07_cell_clustering/", cell_group, "/")
dir.create(results_dir, showWarnings = F, recursive = T)
dir.create(data_dir, showWarnings = F, recursive = T)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

pca_features <- c("TMEM119", "Olig-2", "CD3e", "CD4",
                  "NeuN", "GFAP", "Vimentin", "CD11c",
                  "Iba-1", "CD163", "CD8")
num_pcs <- length(pca_features) - 1
leiden_res <- c(0.1, 0.3, 0.5, 1)
z_cap <- 10
snn_k <- 15
# clustering input: "section_z" = SPACEc (raw, z per section over all cells in it, no cap);
# "log_section_z_group" = PR 32 (log1p, z per section within this group, capped at z_cap)
clust_norm <- "section_z"
stopifnot(clust_norm %in% c("section_z", "log_section_z_group"))

sce <- readRDS("data/05_object_creation/spe.rds")

# save filtered table for now
if (cell_group == names(cell_groups)[1]) {
  arrow::write_parquet(as.data.frame(colData(sce[, sce$keep %in% TRUE]))[, c("code", "cell_id", "cell_type", "area_um2", "region", "has_dpr",
                                                                             "has_ptdp", "dist_dpr_um", "dist_ptdp_um", "dist_vessel_um",
                                                                             "sample_id", "tissue", "sex", "age", "clinical_diagnosis",
                                                                             "c9orf72_mutation", "group", "batch")],
                       paste0("data/05_object_creation/filtered_cells.parquet"))
}

missing <- setdiff(pca_features, rownames(sce))
if (length(missing)) stop("pca_features not in rownames(sce): ", paste(missing, collapse = ", "))

if (clust_norm == "section_z") {
  by_code <- split(seq_len(ncol(sce)), sce$code)
  cnt <- assay(sce, "counts")[pca_features, , drop = F]
  sec_mu <- sapply(by_code, function(i) rowMeans(cnt[, i, drop = F]))
  sec_sd <- sapply(by_code, function(i) matrixStats::rowSds(cnt[, i, drop = F]))
  sec_sd[is.na(sec_sd) | sec_sd == 0] <- 1
  write.csv(rbind(cbind(stat = "mean", marker = pca_features, as.data.frame(sec_mu)),
                  cbind(stat = "sd", marker = pca_features, as.data.frame(sec_sd))),
            paste0(results_dir, "section_z_params.csv"), row.names = F)
  rm(cnt, by_code)
}

missing_types <- setdiff(cell_groups[[cell_group]], unique(sce$cell_type))
if (length(missing_types)) stop("cell types not in the SPE: ", paste(missing_types, collapse = ", "))
uncovered <- setdiff(unique(sce$cell_type), unlist(cell_groups))
if (length(uncovered)) warning("cell types in no group (not embedded): ", paste(uncovered, collapse = ", "))

sce <- sce[, sce$cell_type %in% cell_groups[[cell_group]]]
gc()
assays(sce) <- assays(sce)[c("counts", "asinh")]
sce$cell_type <- factor(sce$cell_type)
sce$code <- factor(sce$code)

clust <- assay(sce, "counts")[pca_features, , drop = F]
if (clust_norm == "section_z") {
  cd <- as.character(sce$code)
  clust <- (clust - sec_mu[, cd, drop = F]) / sec_sd[, cd, drop = F]
} else {
  clust <- log1p(pmax(clust, 0))
  for (cd in levels(sce$code)) {
    i <- which(sce$code == cd)
    mu <- rowMeans(clust[, i, drop = F])
    sdv <- matrixStats::rowSds(clust[, i, drop = F])
    sdv[is.na(sdv) | sdv == 0] <- 1
    clust[, i] <- pmax(pmin((clust[, i, drop = F] - mu) / sdv, z_cap), -z_cap)
  }
}
assay(sce, "counts") <- NULL
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

reducedDim(sce, "PCA") <- scater::calculatePCA(clust, ncomponents = length(pca_features),
                                               ntop = nrow(clust), scale = F)
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

# Leiden clustering -------------------------------------------------------

snn_args <- list(reducedDim(sce, "HARMONY")[, 1:num_pcs], k = snn_k, BNPARAM = BiocNeighbors::AnnoyParam())
if ("num.threads" %in% names(formals(bluster::makeSNNGraph))) snn_args$num.threads <- n_cores else
  snn_args$BPPARAM <- BiocParallel::MulticoreParam(n_cores)
set.seed(1)
graph <- do.call(bluster::makeSNNGraph, snn_args)

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
    if (v == "cell_type" && nlevels(sce$cell_type) < 2) next
    ct <- table(cluster = cl, colData(sce)[[v]])
    write.csv(as.data.frame.matrix(ct), paste0(results_dir, col, "_", v, "_counts.csv"))
    frac_heatmap(ct, paste0(col, "_", v, "_fraction.png"), paste("Fraction of each cluster by", v))
  }

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
                    main = "AUROC, cluster vs rest of group", silent = T)
  h_lfc <- pheatmap(pmax(pmin(lfc_m[ord, ], lim), -lim), cluster_rows = F, cluster_cols = F, color = pal,
                    breaks = seq(-lim, lim, length.out = 102), fontsize_row = 7,
                    main = sprintf("log2FC vs rest (mean asinh diff / ln 2; capped at %.2f)", lim), silent = T)
  png(paste0(results_dir, col, "_markers_auc_log2fc.png"), units = "in", res = 300,
      height = 12, width = 2 * max(5, nlevels(cl) * 0.3 + 3), bg = "white")
  gridExtra::grid.arrange(h_auc$gtable, h_lfc$gtable, ncol = 2)
  dev.off()
}
