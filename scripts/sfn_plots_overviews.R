dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
suppressPackageStartupMessages({
  library(tidyverse)
  library(SpatialExperiment)
  library(scater)
  library(presto)
  library(dittoSeq)
  library(harmony)
  library(pheatmap)
  library(speckle)
})

setwd("/projects/b1169/boles/als_cns_pcf")

dims <- readRDS("data/09_qc_dimred/dataset_z_all_features/reduced_dims.rds")
umap <- dims$HARMONY_UMAP   # raw UMAP (no Harmony); use dims$HARMONY_UMAP for the Harmony one

meta <- read_csv("data/09_qc_dimred/filtered_data.csv", col_select = c(cell_id, cell_type2, code))

out <- tibble(cell_id = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2]) %>%
  left_join(meta, by = "cell_id")
stopifnot(!anyNA(out$cell_type2))

write_csv(out, "results/09_qc_dimred/asinh/raw_umap_coords.csv")

# UMAP --------------------------------------------------------------------

out <- read.csv("/projects/b1169/boles/als_cns_pcf/results/09_qc_dimred/asinh/raw_umap_coords.csv")

source("/projects/b1169/boles/als_cns_scrnaseq/r_scripts/plot_theme.R")

fig1_pal <- DiscretePalette_scCustomize(24, "polychrome")[c(2, 7, 10, 13, 16, 24)]

celltype_cols <- c(
  Astrocyte       = "#39866C",  # pale gold
  Macrophage      = "#8E3B2E",  # dark brick
  Microglia       = "#E07A2E",  # orange
  Neuron          = "#4C439A",  # indigo (unchanged)
  Oligodendrocyte = "#7CAEDB",  # light blue (unchanged)
  `T-cell`        = "#E5B975"   # rose (unchanged)
)

p <- out %>% 
  # slice_sample(n = 20000, by = cell_type2) %>%
  slice_sample(prop = 1) %>%
  mutate(cell_type2 = factor(cell_type2,
                             levels = c("astrocyte", "Macrophage", "Microglia", "cortical_neuron", "oligodendrocyte", "cd4_tcell", "cd8_tcell"),
                             labels = c("Astrocyte", "Macrophage", "Microglia", "Neuron", "Oligodendrocyte", "T-cell", "T-cell"))) %>%
  ggplot(aes(x = UMAP_1,
             y = UMAP_2)) + 
  ggrastr::rasterize(
  geom_point(aes(color = cell_type2),
             shape = 16,
             size = 0.5),
  dpi = 1600
  ) + 
  scale_color_manual(values = celltype_cols) +
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 2))) + 
  umap_labs + 
  theme_classic() +
  theme_als_project(base_size = 8) +
  umap_theme + 
  theme(plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(margin = margin(r = 0.01, unit = "cm")),
        # legend.key.spacing.y = unit(0.01, "cm"),
        legend.box.spacing = unit(0.1, "cm"))
ggsave(p,
       filename = "results/sfn_figs/umap.png",
       units = "in", dpi = 600,
       height = 4, width = 5)


# Bar chart? --------------------------------------------------------------

df <- read.csv("data/09_qc_dimred/filtered_data.csv")

df <- df %>%
  mutate(cell_type2 = factor(cell_type2,
                             levels = c("astrocyte", "Macrophage", "Microglia", "cortical_neuron", "oligodendrocyte", "cd4_tcell", "cd8_tcell"),
                             labels = c("Astrocyte", "Macrophage", "Microglia", "Neuron", "Oligodendrocyte", "T-cell", "T-cell")))

str(df)

cell_totals <- df %>% 
  filter(region %in% c("GM", "WM")) %>%
  group_by(region, group) %>%
  summarize(total = n())

cell_counts <- df %>% 
  filter(region %in% c("GM", "WM")) %>% 
  group_by(region, group, cell_type2) %>% 
  summarize(count = n())

full_join(cell_counts, cell_totals,
          by = c("region", "group")) %>% 
  mutate(prop = count / total) %>% 
  mutate(group = factor(group, levels = c("Control", "sALS", "C9orf72-ALS")),
         region = factor(region,
                         levels = c("GM", "WM"),
                         labels = c("Gray matter", "White matter"))) %>%
  ggplot(aes(x = group,
             y = prop)) + 
  geom_col(aes(fill = cell_type2),
           color = "black") + 
  scale_fill_manual(values = celltype_cols) +
  facet_wrap(. ~ region,
             ncol = 1) + 
  labs(y = "Proportion of all cells") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_linedraw() +
  theme_als_project(base_size = 8) + 
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
ggsave(last_plot(),
       filename = "results/sfn_figs/cell_props_bars.png",
       units = "in", dpi = 600,
       height = 4, width = 3)

propeller()
# Spatial plot? -----------------------------------------------------------

str(df)

p <- df %>% 
  filter(code == "162-1") %>% 
  slice_sample(prop = 1) %>%
  ggplot(aes(x = x,
             y = y)) + 
  ggrastr::rasterize(
    geom_point(aes(color = cell_type2),
               shape = 16,
               size = 0.1),
    dpi = 1600
  ) + 
  scale_color_manual(values = celltype_cols) +
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4))) + 
  theme_als_project(base_size = 8) +
  theme_void() +
  theme(plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(margin = margin(r = 0.01, unit = "cm")),
        # legend.key.spacing.y = unit(0.01, "cm"),
        legend.box.spacing = unit(0.1, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        aspect.ratio = 1)
ggsave(p,
       filename = "results/sfn_figs/162-1_spatial_allcells.png",
       units = "in", dpi = 600,
       height = 4, width = 5)
