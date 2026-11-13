dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
suppressPackageStartupMessages({
  library(tidyverse)
  library(SpatialExperiment)
  library(scater)
  library(dittoSeq)
  library(paletteer)
  library(speckle)
})
setwd("/projects/b1169/boles/als_cns_pcf")

source("/projects/b1169/boles/als_cns_scrnaseq/r_scripts/plot_theme.R")

run <- "log_section_z"

dims <- readRDS(paste0("data/11_myeloid_dimred/", run, "/reduced_dims.rds"))
clusters <- arrow::read_parquet(paste0("data/12_myeloid_clustering/", run, "/myeloid/clusters.parquet"))
labels <- data.table::fread("data/09_qc_dimred/filtered_data.csv", select = c("cell_id", "cell_type2"), data.table = F)

spe <- readRDS("data/05_object_creation/spe.rds")
spe <- spe[, match(rownames(dims$PCA), spe$cell_id)]           # 11's cells, in 11's order
stopifnot(identical(spe$cell_id, rownames(dims$PCA)))

for (d in names(dims)) reducedDim(spe, d) <- dims[[d]]         # PCA, UMAP, SUP_UMAP, HARMONY, HARMONY_UMAP, HARMONY_SUP_UMAP

spe$cell_type2 <- factor(labels$cell_type2[match(spe$cell_id, labels$cell_id)])
cl <- clusters[match(spe$cell_id, clusters$cell_id), ]
stopifnot(identical(cl$cell_id, spe$cell_id))
for (col in grep("^leiden_", names(cl), value = T)) colData(spe)[[col]] <- factor(cl[[col]])

rm(dims, clusters, labels, cl); gc()
spe

# quick checks
table(spe$cell_type2, spe$leiden_0.3)
dittoDimPlot(spe, "leiden_0.3", reduction.use = "HARMONY_UMAP", do.label = T, order = "randomize")

# UMAP of Leiden res = 0.3 ------------------------------------------------

dittoDimPlot(spe,
             "leiden_0.3",
             reduction.use = "HARMONY_UMAP",
             do.label = F,
             order = "randomize",
             color.panel = paletteer_d("ggsci::default_locuszoom")) + 
  umap_labs + 
  theme_als_project(base_size = 8) + 
  umap_theme +
  theme(plot.title = element_blank(),
        legend.box.spacing = unit(0.1, "cm"))
ggsave(last_plot(),
       filename = "results/sfn_figs/myeloid_umap.png",
       units = "in", dpi = 600,
       height = 4, width = 4.5)

spe$group <- factor(spe$group,
                   levels = c("Control", "sALS", "C9orf72-ALS"))

dittoDimPlot(spe,
             "group",
             reduction.use = "HARMONY_UMAP",
             do.label = F,
             order = "randomize") + 
  umap_labs + 
  theme_als_project(base_size = 8) + 
  umap_theme +
  theme(plot.title = element_blank(),
        legend.box.spacing = unit(0.1, "cm"))

# Cluster changes by group ------------------------------------------------

idx <- spe$region %in% c("GM", "WM")

propeller(clusters = spe$leiden_0.3[idx],
          sample = spe$code[idx],
          group = spe$group[idx])

dittoBarPlot(spe,
             var = "group",
             group.by = "leiden_0.5")

dittoPlot(spe,
          var = "GPNMB",
          group.by = "group",
          # plots = "violin",
          assay = "asinh")

dittoDimPlot(spe,
             "region",
             reduction.use = "HARMONY_UMAP")

microglia_ids <- colnames(spe[, spe$cell_type2 == "Microglia"])

de_c9 <- wilcoxauc(spe, group_by = "group",
          assay = "asinh",
          cells.use = microglia_ids,
          groups_use = c("Control", "C9orf72-ALS"))

de_sals <- wilcoxauc(spe, group_by = "group",
                   assay = "asinh",
                   cells.use = microglia_ids,
                   groups_use = c("Control", "sALS"))


# de_c9 %>% 
#   filter(group == "C9orf72-ALS") %>%
#   # mutate(log2FC = logFC / log(2)) %>%
#   arrange(desc(auc)) %>% 
#   slice_head(n = 10)

de_c9 <- de_c9 %>% 
  filter(group != "Control") %>% 
  dplyr::select(c(feature, logFC, auc)) %>% 
  mutate(logFC = logFC / log(2)) %>%
  dplyr::rename("c9_fc" = "logFC", "c9_auc" = "auc")

de_sals <- de_sals %>% 
  filter(group != "Control") %>% 
  dplyr::select(c(feature, logFC, auc)) %>% 
  mutate(logFC = logFC / log(2)) %>%
  dplyr::rename("sals_fc" = "logFC", "sals_auc" = "auc")

# de_sals %>% 
#   filter(group == "sALS") %>%
#   arrange(desc(auc)) %>% 
#   slice_head(n = 10)

de_full <- full_join(de_c9,
                     de_sals,
                     by = "feature")

label_features <- de_full %>%
  mutate(magnitude = sqrt(sals_fc^2 + c9_fc^2),
         sig_group = if_else(c9_fc > 0 | sals_fc > 0, "up", "down")) %>%
  group_by(sig_group) %>%
  slice_max(magnitude, n = 5) %>%
  ungroup()

de_full %>% 
  ggplot(aes(x = sals_fc,
             y = c9_fc)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray60") + 
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "gray60") +
  geom_point(fill = "dodgerblue",
             shape = 21,
             size = 2) + 
  geom_label_repel(data = label_features,
                   aes(label = feature),
                   show.legend = FALSE,
                   seed = 68,                    # reproducible layout
                   force = 25,
                   force_pull = 0.5,
                   max.overlaps = Inf,            # never drop a label
                   box.padding = 0,             # space around label
                   point.padding = 0,           # space around the point it's labeling
                   min.segment.length = 0,        # always draw the leader line
                   segment.color = "grey40",
                   # THIS is the part that avoids all points, not just labeled ones:
                   max.time = 5,
                   max.iter = 100000,
                   size = 2.5) + 
  ggtitle("Microglial protein markers in ALS") +
  labs(y = "FC in C9orf72-ALS",
       x = "FC in sALS") +
  theme_linedraw() + 
  theme_als_project(base_size = 8)
ggsave(last_plot(),
       filename = "results/sfn_figs/myeloid_markers.png",
       units = "in", dpi = 600,
       height = 6, width = 6)

dittoPlot(spe[, spe$group != "Control" & spe$region == "GM" & spe$cell_type2 == "Microglia"],
          "dist_ptdp_um",
          # plots = "ridgeplot",
          plots = c("vlnplot", "jitter"),
          group.by = "group",
          color.panel = c("#0CAA00", "#CC00FF")) + 
  ggtitle("Distance to nearest\nTDP-43 inclusion") +
  labs(y = "Microns") +
  theme_als_project(base_size = 8) + 
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave(last_plot(),
       filename = "results/sfn_figs/myeloid_tdp_dist.png",
       units = "in", dpi = 600,
       height = 4, width = 3)


dist_tab <- tibble(
  id = spe$cell_id,
  type = spe$cell_type2,
  region = spe$region,
  group = spe$group,
  tdp = spe$dist_ptdp_um,
  dpr = spe$dist_dpr_um,
) %>% 
  filter(region == "GM" &
           group == "C9orf72-ALS" & 
           type == "Microglia")

dist_tab %>% 
  pivot_longer(c(tdp, dpr),
               names_to = "feature",
               values_to = "distance") %>% 
  ggplot(aes(x = feature,
             y = distance)) + 
  geom_violin() + 
