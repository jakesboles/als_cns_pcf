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

# Get myeloid subtypes from 09 --------------------------------------------

