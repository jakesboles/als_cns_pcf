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

# transform, then per-marker z over the whole dataset or per section (code); one array task per method
methods <- list(
  asinh           = list(transform = "asinh", scope = "dataset", cap = Inf, harmony = T),
  dataset_z       = list(transform = "none",  scope = "dataset", cap = Inf, harmony = T),
  section_z       = list(transform = "none",  scope = "section", cap = Inf, harmony = T),
  log_dataset_z   = list(transform = "log1p", scope = "dataset", cap = 10,  harmony = T),
  log_section_z   = list(transform = "log1p", scope = "section", cap = 10,  harmony = T),
  asinh_section_z = list(transform = "asinh", scope = "section", cap = 10,  harmony = T),
  # AD CODEX-CNS: log1p -> ComBat by section -> sc.pp.scale(max_value = 10); no Harmony
  ad_combat       = list(transform = "log1p", scope = "combat",  cap = 10,  harmony = F)
)

method <- Sys.getenv("PCF_METHOD")
if (method == "") {
  task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
  if (is.na(task) || task < 1 || task > length(methods))
    stop("Set PCF_METHOD or run as array task 1-", length(methods), ": ", paste(names(methods), collapse = ", "))
  method <- names(methods)[task]
}
if (!method %in% names(methods)) stop("Unknown PCF_METHOD ", method, "; one of ", paste(names(methods), collapse = ", "))

# panel: pca_features; all: every channel kept in 05
feature_set <- Sys.getenv("PCF_FEATURES", "panel")
if (!feature_set %in% c("panel", "all")) stop("PCF_FEATURES must be panel or all")
run_name <- paste0(method, if (feature_set == "all") "_all_features" else "")
message("Method: ", method, "; features: ", feature_set, "; outputs: ", run_name)

results_dir <- paste0("results/11_myeloid_dimred/", run_name, "/")
data_dir <- paste0("data/11_myeloid_dimred/", run_name, "/")
dir.create(results_dir, showWarnings = F, recursive = T)
dir.create(data_dir, showWarnings = F, recursive = T)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# microglia / macrophage panel
pca_features <- c("Iba-1", "TMEM119", "CD163", "CD11c", "HLA-DR", "CD68", "CD74", "GPNMB", "Galectin-3",
                  "LYVE1", "SPP1", "STING", "MMP9", "Vimentin", "APOC1", "ApoE", "PSAP", "CD14")
npcs_env <- Sys.getenv("PCF_NPCS")
# supervised UMAP (display only): weight of cell_type2 vs the data, 0-1
sup_weight <- as.numeric(Sys.getenv("PCF_SUP_WEIGHT", "0.5"))

sce <- readRDS("data/05_object_creation/spe.rds")

# Myeloid cells kept by 09 -------------------------------------------------

keep <- data.table::fread("data/09_qc_dimred/filtered_data.csv", select = c("cell_id", "cell_type2"), data.table = F) %>%
  filter(cell_type2 %in% c("Microglia", "Macrophage"))
sce <- sce[, sce$cell_id %in% keep$cell_id]
stopifnot(ncol(sce) == nrow(keep))
sce$cell_type2 <- keep$cell_type2[match(sce$cell_id, keep$cell_id)]
rm(keep)
gc()
print(table(sce$code, sce$cell_type2))

# Pre-processing ----------------------------------------------------------

sce$cell_type2 <- factor(sce$cell_type2)
sce$code <- factor(sce$code)

clust_features <- if (feature_set == "all") rownames(sce) else pca_features
missing <- setdiff(clust_features, rownames(sce))
if (length(missing)) stop("pca_features not in rownames(sce): ", paste(missing, collapse = ", "))
num_pcs <- if (nzchar(npcs_env)) as.integer(npcs_env) else if (feature_set == "all") 12 else length(clust_features) - 1

m <- methods[[method]]
x <- as.matrix(assay(sce, if (m$transform == "asinh") "asinh" else "counts")[clust_features, , drop = F])
if (m$transform == "log1p") x <- log1p(pmax(x, 0))
zscore <- function(v, cap = m$cap, lower = T) {
  mu <- rowMeans(v)
  sdv <- matrixStats::rowSds(v)
  sdv[is.na(sdv) | sdv == 0] <- 1
  v <- pmin((v - mu) / sdv, cap)
  if (lower) pmax(v, -cap) else v
}
ad_display <- NULL
if (m$scope == "dataset") {
  x <- zscore(x)
} else if (m$scope == "combat") {
  if (!requireNamespace("sva", quietly = T)) stop("ad_combat needs the Bioconductor package sva")
  x <- sva::ComBat(x, batch = as.character(sce$code), par.prior = T, prior.plots = F)
  # scanpy's scale clips only the upper end
  x <- zscore(x, lower = F)
  reducedDim(sce, "SCALED") <- t(x)
  # their display UMAP input: z within sample, then z over all cells
  u <- x
  for (cd in levels(sce$code)) {
    i <- which(sce$code == cd)
    u[, i] <- zscore(u[, i, drop = F], cap = Inf)
  }
  ad_display <- t(zscore(u, cap = Inf))
  rm(u)
} else {
  for (cd in levels(sce$code)) {
    i <- which(sce$code == cd)
    x[, i] <- zscore(x[, i, drop = F])
  }
}
assays(sce) <- assays(sce)["asinh"]
gc()
message(length(clust_features), " PCA features, ", num_pcs, " PCs used, ", ncol(sce), " cells")

# Dimensional reduction ---------------------------------------------------

reducedDim(sce, "PCA") <- scater::calculatePCA(x, scale = F, ncomponents = min(length(clust_features), 30))
rm(x)
gc()
stopifnot(num_pcs <= ncol(reducedDim(sce, "PCA")))

rotation <- attr(reducedDim(sce, "PCA"), "rotation")
percent_var <- attr(reducedDim(sce, "PCA"), "percentVar")
write.csv(rotation, paste0(results_dir, "pca_loadings.csv"))

png(paste0(results_dir, "pca_elbow.png"), units = "in", res = 600, height = 4, width = 5, bg = "white")
plot(percent_var, log = "y", xlab = "PC", ylab = "Variance explained (%)", main = run_name)
dev.off()

png(paste0(results_dir, "pca_loadings_heatmap.png"), units = "in", res = 300, height = max(14, nrow(rotation) * 0.2 + 2),
    width = max(7, ncol(rotation) * 0.35 + 3), bg = "white")
pheatmap(rotation,
         color = colorRampPalette(c("dodgerblue", "white", "firebrick"))(101),
         breaks = seq(-max(abs(rotation)), max(abs(rotation)), length.out = 102),
         cluster_cols = F,
         labels_col = paste0(colnames(rotation), " (", round(percent_var, 1), "%)"),
         fontsize_row = 7, main = run_name)
dev.off()

set.seed(1)
explain_idx <- sample(ncol(sce), min(ncol(sce), 2e5))
explain_vars <- c("code", "batch", "cell_type2", "group", "sex", "age", "area_um2", "region")
explain_vars <- explain_vars[map_lgl(explain_vars, ~ length(unique(na.omit(colData(sce)[[.x]][explain_idx]))) > 1)]

plot_explanatory <- function(dimred, file) {
  ex <- getExplanatoryPCs(sce[, explain_idx], variables = explain_vars, dimred = dimred)
  png(paste0(results_dir, file), units = "in", res = 600, height = 6, width = 7, bg = "white")
  print(plotExplanatoryPCs(ex / 100))
  dev.off()
}

umap_plots <- function(dimred, prefix) {
  files <- c(cell_type2 = "celltype", batch = "batch", code = "sample", group = "group", region = "region")
  for (v in names(files)) {
    p <- dittoDimPlot(sce, v, reduction.use = dimred, order = "randomize", main = paste(run_name, prefix))
    ggsave(paste0(results_dir, prefix, "_umap_", files[[v]], ".png"), p, units = "in", dpi = 600,
           height = 5, width = if (v == "code") 9 else 6)
  }
  marker_plots <- map(pca_features,
                      ~ dittoDimPlot(sce[, explain_idx], .x, assay = "asinh", reduction.use = dimred,
                                     order = "increasing", size = 0.3, legend.show = F) +
                        labs(title = .x, x = NULL, y = NULL))
  p <- gridExtra::arrangeGrob(grobs = marker_plots, ncol = 4)
  ggsave(paste0(results_dir, prefix, "_umap_markers.png"), p, width = 16,
         height = 4 * ceiling(length(marker_plots) / 4), dpi = 300)
}

plot_explanatory("PCA", "explanatory_pcs.png")

set.seed(1)
reducedDim(sce, "UMAP") <- uwot::umap(reducedDim(sce, "PCA")[, 1:num_pcs], n_neighbors = 15,
                                      n_threads = n_cores, n_sgd_threads = n_cores, verbose = T)
umap_plots("UMAP", "raw")

supervised_umap <- function(dimred) {
  set.seed(1)
  uwot::umap(reducedDim(sce, dimred)[, 1:num_pcs], n_neighbors = 15, y = sce$cell_type2,
             target_weight = sup_weight, n_threads = n_cores, n_sgd_threads = n_cores)
}
reducedDim(sce, "SUP_UMAP") <- supervised_umap("PCA")
umap_plots("SUP_UMAP", "raw_supervised")

if (m$harmony) {
  sce <- RunHarmony(sce, group.by.vars = "code")
  plot_explanatory("HARMONY", "explanatory_harmony.png")

  set.seed(1)
  reducedDim(sce, "HARMONY_UMAP") <- uwot::umap(reducedDim(sce, "HARMONY")[, 1:num_pcs], n_neighbors = 15,
                                                n_threads = n_cores, n_sgd_threads = n_cores)
  umap_plots("HARMONY_UMAP", "harmony")
  reducedDim(sce, "HARMONY_SUP_UMAP") <- supervised_umap("HARMONY")
  umap_plots("HARMONY_SUP_UMAP", "harmony_supervised")
}

if (!is.null(ad_display)) {
  # umap-learn defaults, as in their fabian_style_clustering()
  set.seed(1)
  reducedDim(sce, "AD_UMAP") <- uwot::umap(ad_display, n_neighbors = 15, min_dist = 0.1,
                                           n_threads = n_cores, n_sgd_threads = n_cores)
  rm(ad_display)
  umap_plots("AD_UMAP", "ad_display")
}

dr_names <- intersect(c("PCA", "UMAP", "SUP_UMAP", "HARMONY", "HARMONY_UMAP", "HARMONY_SUP_UMAP", "SCALED", "AD_UMAP"),
                      reducedDimNames(sce))
dims <- map(set_names(dr_names), ~ `rownames<-`(reducedDim(sce, .x), sce$cell_id))
saveRDS(dims, paste0(data_dir, "reduced_dims.rds"))

# Separation scores -------------------------------------------------------

# silhouette on up to 5000 cells per cell_type2: high for cell_type2 = types apart; near 0 for code = sections mixed
set.seed(1)
sil_idx <- unlist(map(split(seq_len(ncol(sce)), sce$cell_type2), ~ if (length(.x) > 5000) sample(.x, 5000) else .x))
sil <- map_dfr(intersect(c("PCA", "HARMONY", "SCALED"), reducedDimNames(sce)), function(dr) {
  emb <- reducedDim(sce, dr)
  d <- dist(emb[sil_idx, if (dr == "SCALED") seq_len(ncol(emb)) else 1:num_pcs, drop = F])
  map_dfr(c("cell_type2", "code"), function(v) {
    lab <- colData(sce)[[v]][sil_idx]
    s <- cluster::silhouette(as.integer(factor(lab)), d)[, "sil_width"]
    tibble(method = run_name, embedding = dr, label = v, level = c("all", levels(factor(lab))),
           mean_silhouette = c(mean(s), tapply(s, factor(lab), mean)))
  })
})
write.csv(sil, paste0(results_dir, "silhouette.csv"), row.names = F)
print(filter(sil, level == "all"))
