dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
library(tidyverse)
library(SpatialExperiment)

setwd("/projects/b1169/boles/als_cns_pcf")

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

table(df$leiden_1)

# cell_clustering - res 1 - clusters 9/17 - 60262 - worst pending below
# cell_clustering_allcell_z - res 1 - clusters 10/7 - 76073 - 2nd best pending below
# cell_clustering_group_panel_allcell_z - res 1 - clusters 7 - 66033 - best pending below
# cell_clustering_all_features - res 1 - clusters 7 - 33360 - not better than above
# myeloid subtyping from 09 - res 1 - cluster 5 macrophages, cluster 8 most likely monocytes - 59785/5567 respectively
    # this looks bad! bad segmentations in macrophages! 
