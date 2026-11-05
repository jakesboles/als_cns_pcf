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

explain_pcs <- getExplanatoryPCs(sce,
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
sce <- scater::runUMAP(sce,
               dimred = "PCA",
               ntop = num_pcs)

save_umap <- function(file, w){
  ggsave(filename = paste0(results_dir, file),
         plot = last_plot(),
         units = "in", dpi = 600,
         height = 5, width = w)
}

dittoDimPlot(sce,
             "cell_type",
             reduction.use = "UMAP")
save_umap("raw_umap_celltype.png", 6)

dittoDimPlot(sce,
             "batch",
             reduction.use = "UMAP")
save_umap("raw_umap_batch.png", 6)

dittoDimPlot(sce,
             "code",
             reduction.use = "UMAP")
save_umap("raw_umap_sample.png", 9)

# Harmony -----------------------------------------------------------------

# sce <- RunHarmony(sce,
#                   group.by.vars = "code")
# 
# explain_pcs <- getExplanatoryPCs(sce,
#                                  variables = c("code",
#                                                "batch",
#                                                "group",
#                                                "sex", "age",
#                                                "area"),
#                                  dimred = "HARMONY"
# )
# 
# plotExplanatoryPCs(explain_pcs/100)
# 
# sce <- scater::runUMAP(sce,
#                        dimred = "HARMONY",
#                        ntop = num_pcs,
#                        name = "HARMONY_UMAP")
# 
# sce
# 
# scater::plotUMAP(sce,
#                  colour_by = "code",
#                  dimred = "HARMONY_UMAP")