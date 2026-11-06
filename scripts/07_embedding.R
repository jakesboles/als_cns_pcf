dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
library(tidyverse)
library(SpatialExperiment)
library(scater)
library(presto)
library(dittoSeq)
library(harmony)
library(pheatmap)

setwd("/projects/b1169/boles/als_cns_pcf")

results_dir <- "results/07_embedding/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)
data_dir <- "data/07_embedding/"
dir.create(data_dir, showWarnings = F, recursive = T)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

sce <- readRDS("data/05_object_creation/spe.rds")
sce <- sce[, sce$keep]
sce

# save filtered table for now
arrow::write_parquet(as.data.frame(colData(sce))[, c("code", "cell_id", "cell_type", "area_um2", "region", "has_dpr",
                                                     "has_ptdp", "dist_dpr_um", "dist_ptdp_um", "dist_vessel_um",
                                                     "sample_id", "tissue", "sex", "age", "clinical_diagnosis",
                                                     "c9orf72_mutation", "group", "batch")],
                     paste0("data/05_object_creation/filtered_cells.parquet"))

# Dimensional reduction ---------------------------------------------------

sce <- scater::runPCA(sce,
                      scale = T,
                      exprs_values = "asinh")
sce

# scater::plotPCA(sce,
#         colour_by = "group")
# 
# scater::plotPCA(sce,
#         colour_by = "batch")
# 
# scater::plotPCA(sce,
#         colour_by = "code")

percent.var <- attr(reducedDim(sce), "percentVar")

png(filename = paste0(results_dir, "pca_elbow.png"),
    units = "in", res = 600,
    height = 4, width = 5,
    bg = "white")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
dev.off()

num_pcs <- 12

# plotReducedDim(sce, dimred = "PCA", ncomponents = 3, colour_by = "sample_id")

set.seed(1)
explain_idx <- sample(ncol(sce), min(ncol(sce), 2e5))
explain_pcs <- getExplanatoryPCs(sce[, explain_idx],
                                 variables = c("code",
                                               "batch",
                                               "cell_type",
                                               "group",
                                               "sex", "age",
                                               "area_um2")
)

png(filename = paste0(results_dir, "explanatory_pcs.png"),
    units = "in", res = 600,
    height = 6, width = 7,
    bg = "white")
plotExplanatoryPCs(explain_pcs/100)
dev.off()
# code/batch explain the most variance - integration likely needed

# check if integration is needed with UMAP
set.seed(1)
reducedDim(sce, "UMAP") <- uwot::umap(reducedDim(sce, "PCA")[, 1:num_pcs],
                                      n_neighbors = 15,
                                      n_threads = n_cores,
                                      n_sgd_threads = n_cores,
                                      verbose = T)

# saveRDS(list(PCA = reducedDim(sce, "PCA"), UMAP = reducedDim(sce, "UMAP")),
#         paste0(data_dir, "reduced_dims.rds"))

save_umap <- function(p, file, w){
  ggsave(filename = paste0(results_dir, file),
         plot = p,
         units = "in", dpi = 600,
         height = 5, width = w)
}

p <- dittoDimPlot(sce,
                  "cell_type",
                  reduction.use = "UMAP")
save_umap(p, "raw_umap_celltype.png", 6)

p <- dittoDimPlot(sce,
                  "batch",
                  reduction.use = "UMAP")
save_umap(p, "raw_umap_batch.png", 6)

p <- dittoDimPlot(sce,
                  "code",
                  reduction.use = "UMAP")
save_umap(p, "raw_umap_sample.png", 9)
# desperately need integration! 

# Harmony -----------------------------------------------------------------

sce <- RunHarmony(sce,
                  group.by.vars = "code")

explain_pcs <- getExplanatoryPCs(sce,
                                 variables = c("code",
                                               "batch",
                                               "cell_type",
                                               "group",
                                               "sex", "age",
                                               "area_um2"),
                                 dimred = "HARMONY"
)

png(filename = paste0(results_dir, "explanatory_harmony.png"),
    units = "in", res = 600,
    height = 6, width = 7,
    bg = "white")
plotExplanatoryPCs(explain_pcs/100)
dev.off()

reducedDim(sce, "HARMONY_UMAP") <- uwot::umap(reducedDim(sce, "HARMONY")[, 1:num_pcs],
                                              n_neighbors = 15,
                                              n_threads = n_cores,
                                              n_sgd_threads = n_cores)

sce

saveRDS(list(HARMONY = reducedDim(sce, "HARMONY"), 
             HARMONY_UMAP = reducedDim(sce, "HARMONY_UMAP")),
        paste0(data_dir, "reduced_dims_harmony.rds"))

p <- dittoDimPlot(sce,
                  "cell_type",
                  reduction.use = "HARMONY_UMAP")
save_umap(p, "harmony_umap_celltype.png", 6)

p <- dittoDimPlot(sce,
                  "batch",
                  reduction.use = "HARMONY_UMAP")
save_umap(p, "harmony_umap_batch.png", 6)

p <- dittoDimPlot(sce,
                  "code",
                  reduction.use = "HARMONY_UMAP")
save_umap(p, "harmony_umap_sample.png", 9)