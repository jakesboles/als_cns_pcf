dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
library(tidyverse)
library(SpatialExperiment)
library(scater)
library(presto)
library(dittoSeq)
library(harmony)
library(pheatmap)

setwd("/projects/b1169/boles/als_cns_pcf")

sce <- readRDS("data/05_object_creation/spe.rds")
sce <- sce[, sce$keep]
sce

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
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
num_pcs <- 10

# plotReducedDim(sce, dimred = "PCA", ncomponents = 3, colour_by = "sample_id")

explain_pcs <- getExplanatoryPCs(sce,
                                 variables = c("code",
                                               "batch",
                                               "cell_type",
                                               # "group",
                                               # "sex", "age",
                                               "area_um2")
)

plotExplanatoryPCs(explain_pcs/100)

sce <- scater::runUMAP(sce,
               dimred = "PCA",
               ntop = num_pcs)

scater::plotUMAP(sce,
         colour_by = "cell_type")

dittoDimPlot(sce,
             "NeuN",
             assay = "asinh",
             reduction.use = "UMAP")

dittoDimPlot(sce,
             "dist_vessel_um",
             reduction.use = "UMAP")

dittoDimPlot(sce,
             "has_dpr",
             reduction.use = "UMAP")

sce$touching_vessel <- if_else(sce$dist_vessel_um == 0, "Y", "N")

dittoDimPlot(sce,
             "touching_vessel",
             reduction.use = "UMAP")

# 
# scater::plotUMAP(sce,
#          colour_by = "batch")

# scater::plotUMAP(sce,
#          colour_by = "code")

# Harmony -----------------------------------------------------------------

sce <- RunHarmony(sce,
                  group.by.vars = "code")

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

sce <- scater::runUMAP(sce,
                       dimred = "HARMONY",
                       ntop = num_pcs,
                       name = "HARMONY_UMAP")

sce

scater::plotUMAP(sce,
                 colour_by = "code",
                 dimred = "HARMONY_UMAP")