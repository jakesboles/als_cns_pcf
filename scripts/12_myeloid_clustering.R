dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
suppressPackageStartupMessages({
  library(tidyverse)
  library(SpatialExperiment)
  library(scater)
  library(presto)
  library(dittoSeq)
  library(pheatmap)
})

setwd("/projects/b1169/boles/als_cns_pcf")

# 11 output folders to cluster; one array task per entry
runs <- c("log_section_z", "dataset_z")

run_name <- Sys.getenv("PCF_METHOD")
if (run_name == "") {
  task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
  if (is.na(task) || task < 1 || task > length(runs))
    stop("Set PCF_METHOD or run as array task 1-", length(runs), ": ", paste(runs, collapse = ", "))
  run_name <- runs[task]
}
dims_file <- paste0("data/11_myeloid_dimred/", run_name, "/reduced_dims.rds")
if (!file.exists(dims_file)) stop(dims_file, " not found; run 11 for ", run_name)
message("Clustering 11's Harmony embedding for ", run_name)

results_dir <- paste0("results/12_myeloid_clustering/", run_name, "/")
data_dir <- paste0("data/12_myeloid_clustering/", run_name, "/myeloid/")
dir.create(results_dir, showWarnings = F, recursive = T)
dir.create(data_dir, showWarnings = F, recursive = T)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
leiden_res <- c(0.1, 0.2, 0.3, 0.5, 1)
snn_k <- 15

dims <- readRDS(dims_file)
if (!all(c("HARMONY", "HARMONY_UMAP") %in% names(dims))) stop(run_name, " has no Harmony embedding (e.g. ad_combat)")
npcs_env <- Sys.getenv("PCF_NPCS")
# same default as 11: panel features - 1, 12 with all features
num_pcs <- if (nzchar(npcs_env)) as.integer(npcs_env) else if (grepl("_all_features$", run_name)) 12 else ncol(dims$HARMONY) - 1
message(num_pcs, " Harmony dimensions used")

sce <- readRDS("data/05_object_creation/spe.rds")
sce <- sce[, match(rownames(dims$HARMONY), sce$cell_id)]
stopifnot(identical(sce$cell_id, rownames(dims$HARMONY)))
assays(sce) <- assays(sce)["asinh"]
gc()
keep <- data.table::fread("data/09_qc_dimred/filtered_data.csv", select = c("cell_id", "cell_type2"), data.table = F)
sce$cell_type2 <- factor(keep$cell_type2[match(sce$cell_id, keep$cell_id)])
rm(keep)
stopifnot(!anyNA(sce$cell_type2))
sce$code <- factor(sce$code)
reducedDim(sce, "HARMONY") <- dims$HARMONY
reducedDim(sce, "HARMONY_UMAP") <- dims$HARMONY_UMAP
rm(dims)
gc()
print(table(sce$code, sce$cell_type2))

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

# 08-compatible: cell_type column holds cell_type2
arrow::write_parquet(bind_cols(tibble(cell_id = sce$cell_id, cell_type = as.character(sce$cell_type2)),
                               as.data.frame(colData(sce)[, cluster_cols])),
                     paste0(data_dir, "clusters.parquet"))

# Cluster QC --------------------------------------------------------------

frac_heatmap <- function(tab, file, main) {
  png(paste0(results_dir, file), units = "in", res = 300,
      height = max(4, nrow(tab) * 0.2 + 1.5), width = max(6, ncol(tab) * 0.3 + 2), bg = "white")
  pheatmap(prop.table(tab, 1), cluster_rows = F, cluster_cols = F, display_numbers = T,
           number_format = "%.2f", fontsize_number = 5, main = main)
  dev.off()
}

meta_vars <- c(code = "code", batch = "batch", group = "group", region = "region", cell_type2 = "cell_type2")

for (col in cluster_cols) {
  cl <- colData(sce)[[col]]

  p <- dittoDimPlot(sce, col, reduction.use = "HARMONY_UMAP", do.label = T, labels.size = 3,
                    legend.show = F, order = "randomize", main = paste(run_name, col))
  ggsave(paste0(results_dir, col, "_umap.png"), p, width = 7, height = 6, dpi = 300)

  for (v in names(meta_vars)) {
    ct <- table(cluster = cl, colData(sce)[[meta_vars[[v]]]])
    write.csv(as.data.frame.matrix(ct), paste0(results_dir, col, "_", v, "_counts.csv"))
    frac_heatmap(ct, paste0(col, "_", v, "_fraction.png"), paste("Fraction of each cluster by", v))
  }
  frac_heatmap(t(table(cluster = cl, sce$code)), paste0(col, "_code_by_cluster_fraction.png"),
               "Fraction of each section's cells in each cluster")

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
                    main = "AUROC, cluster vs other myeloid cells", silent = T)
  h_lfc <- pheatmap(pmax(pmin(lfc_m[ord, ], lim), -lim), cluster_rows = F, cluster_cols = F, color = pal,
                    breaks = seq(-lim, lim, length.out = 102), fontsize_row = 7,
                    main = sprintf("log2FC vs rest (mean asinh diff / ln 2; capped at %.2f)", lim), silent = T)
  png(paste0(results_dir, col, "_markers_auc_log2fc.png"), units = "in", res = 300,
      height = 12, width = 2 * max(5, nlevels(cl) * 0.3 + 3), bg = "white")
  gridExtra::grid.arrange(h_auc$gtable, h_lfc$gtable, ncol = 2)
  dev.off()
}
