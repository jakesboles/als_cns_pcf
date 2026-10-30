dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
library(tidyverse)
library(ggplot2)
library(janitor)
library(ggbeeswarm)
library(ggpubr)
library(SpatialExperiment)
library(scater)
library(presto)
library(dittoSeq)
library(harmony)
library(pheatmap)

setwd("/projects/b1169/boles/als_cns_pcf")

features <- read.csv("results/01_quantify_cells/160-1_channels.csv")

features <- features %>% 
  pull(channel)

pq1 <- arrow::read_parquet("data/01_quantify_cells/160-1.parquet")
pq2 <- arrow::read_parquet("data/02_annotate_regions/160-1.parquet")

pq <- full_join(pq1, pq2,
                by = c("code", "cell_id"))

# Add in pathology and vasculature relationship data here -----------------

pq3 <- arrow::read_parquet("data/03_annotate_pathology/160-1.parquet")

pq3 <- pq3 %>% 
  dplyr::select(c(cell_id, has_dpr, has_ptdp))

pq4 <- arrow::read_parquet("data/04_feature_distances/160-1.parquet")

pq <- pq %>% 
  left_join(pq3, 
            by = "cell_id") %>% 
  left_join(pq4,
            by = c("code", "cell_id"))

# Add in demographic data -------------------------------------------------

# Make sure Galectin-3 is named the same every time -----------------------

# galectin3_idx <- colnames(pq) %in% "DAPI"

# colnames(pq)[galectin3_idx] <- "Galectin-3"

# table(pq$region)
# 
# pq %>% 
#   filter(region %in% c("GM", "WM")) %>%
#   ggplot(aes(x = x,
#              y = y)) + 
#   geom_point(aes(color = cell_type),
#              size = 0.1) + 
#   facet_wrap(. ~ region,
#              nrow = 1) +
#   theme_void()

# Filter on size ----------------------------------------------------------

pq %>% 
  ggplot(aes(x = area_um2)) + 
  geom_histogram() + 
  facet_wrap(. ~ cell_type,
             ncol = 1,
             scales = "free_y") +
  scale_x_log10()

pq <- pq %>% 
  filter(area_um2 > 40 & 
           area_um2 < 1500)

# Prep SPE object ---------------------------------------------------------
cell_id <- pq$cell_id

meta <- pq %>% 
  dplyr::select(-features) %>%
  mutate(batch = str_split_i(code, "-", i = 1))

counts <- pq %>%
  dplyr::select(features) %>% 
  t()

rownames(meta) <- cell_id -> colnames(counts)

sce <- SpatialExperiment(assays = list(counts = counts),
                         colData = meta,
                         spatialCoordsNames = c("x", "y"))

# dittoHeatmap(sce,
#              genes = features,
#              annot.by = "code",
#              assay = "counts")

# Arcsinh transform -------------------------------------------------------

set.seed(1)
cof_sub <- assay(sce, "counts")[, sample(ncol(sce), min(ncol(sce), 5e4))]

idx <- which(is.na(cof_sub), arr.ind = TRUE)

cbind(
  row_name = rownames(cof_sub)[idx[, 1]],
  col_name = colnames(cof_sub)[idx[, 2]]
)

bg_peak <- function(v) {
  d <- density(log(v[v > 0]), n = 512)
  pk <- which(diff(sign(diff(d$y))) == -2) + 1
  pk <- pk[d$y[pk] > 0.05 * max(d$y)]
  if (!length(pk)) pk <- which.max(d$y)
  exp(d$x[min(pk)])
}

cofactor <- apply(cof_sub, 1, bg_peak)

# manual overrides, e.g. c(Iba1 = 120)
cofactor_manual <- c()
if (length(cofactor_manual)) cofactor[names(cofactor_manual)] <- cofactor_manual
stopifnot(identical(names(cofactor), rownames(sce)))

as.data.frame(t(cof_sub)) %>% 
  pivot_longer(everything(), names_to = "marker") %>% 
  filter(value > 0) %>% 
  ggplot(aes(x = value)) + 
  geom_density() + 
  geom_vline(data = enframe(cofactor, name = "marker", value = "cofactor"),
             aes(xintercept = cofactor),
             color = "firebrick") + 
  scale_x_log10() + 
  facet_wrap(. ~ marker, scales = "free") + 
  theme_linedraw(base_size = 8)

assay(sce, "asinh") <- asinh(assay(sce, "counts") / cofactor)

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