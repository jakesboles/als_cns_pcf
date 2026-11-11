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

results_dir <- paste0("results/09_qc_dimred/", method, "/")
data_dir <- paste0("data/09_qc_dimred/", method, "/")
dir.create(results_dir, showWarnings = F, recursive = T)
dir.create(data_dir, showWarnings = F, recursive = T)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

pca_features <- c("TMEM119", "Olig-2", "CD3e", "CD4",
                  "NeuN", "GFAP", "Vimentin", "CD11c",
                  "Iba-1", "CD163", "CD8")

leiden_res <- c(0.1, 0.3, 0.5, 1)
z_cap <- 10
snn_k <- 15

sce <- readRDS("data/05_object_creation/spe.rds")

# Flag bad clusters from 07 -----------------------------------------------

remove_params <- list(
  astrocyte = list("allcell", 1, 13),
  oligodendrocyte = list("withincell", 1, 11),
  tcell = list("withincell", 0.3, c(2,5))
)

bad_cells <- list()

remove_params

for (i in seq_along(remove_params)){
  cell <- names(remove_params)[i]
  
  if (remove_params[[i]][[1]] == "allcell") {
    in_dir <- paste0("data/07_cell_clustering_allcell_z/", cell, "/")
  } else { 
    in_dir <- paste0("data/07_cell_clustering/", cell, "/")
  }
  
  res <- remove_params[[i]][[2]]
  
  bad_clusters <- remove_params[[i]][[3]]
  
  df <- arrow::read_parquet(paste0(in_dir, "clusters.parquet"))
  
  bad_cells[[i]] <- df %>% 
    filter(!!sym(paste0("leiden_", res)) %in% bad_clusters) %>% 
    pull(cell_id)
}

# Get myeloid subtypes from 07 --------------------------------------------

df <- arrow::read_parquet("data/07_cell_clustering_group_panel_allcell_z/myeloid/clusters.parquet")
# df <- arrow::read_parquet("data/09_myeloid_subtypes/myeloid/clusters.parquet")

df %>% 
  filter(leiden_1 %in% c(7)) %>% nrow()

macrophage_ids <- df %>% 
  filter(leiden_1 %in% c(7)) %>%
  pull(cell_id)

# Update full object and save filtered table for early analysis -----------

meta <- as.data.frame(colData(sce)) %>% 
  mutate(cell_type2 = case_when(cell_id %in% list_c(bad_cells) ~ "Remove",
                                cell_id %in% macrophage_ids ~ "Macrophage",
                                !(cell_id %in% macrophage_ids) & 
                                  cell_type == "myeloid" ~ "Microglia",
                                .default = cell_type))

table(meta$cell_type2, meta$cell_type)

sce$cell_type2 <- meta$cell_type2

table(sce$cell_type2, sce$cell_type)

sce_filtered <- sce[, sce$cell_type2 != "Remove"]

meta <- colData(sce_filtered)[, c("code", "cell_id", "cell_type", "cell_type2", "area_um2", "region", "has_dpr",
                                  "has_ptdp", "dist_dpr_um", "dist_ptdp_um", "dist_vessel_um",
                                  "sample_id", "tissue", "sex", "age", "clinical_diagnosis",
                                  "c9orf72_mutation", "group", "batch")] %>% 
  as.data.frame()
coords <- spatialCoords(sce_filtered)

meta <- cbind(meta, coords)

write.csv(meta,
          file = "data/09_qc_dimred/filtered_data.csv",
          row.names = F)
