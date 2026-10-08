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

setwd("/projects/b1169/boles/als_cns_pcf")

# Prep demographics and joining key ---------------------------------------

key <- read.csv("sample_key.csv")

key <- key[, c("code", "sample_id", "tissue")]

demographics <- read.csv("target_als_demographics_compiled.csv")

demographics <- demographics[, c("Case.Number", "Sex", "Age.at.Death", "Clinical.Diagnosis", "C9orf72.mutation")] %>%
  mutate(group = case_when(C9orf72.mutation == "Y" ~ "C9orf72-ALS",
                           Clinical.Diagnosis == "Control" ~ "Control",
                           Clinical.Diagnosis != "Control" & C9orf72.mutation != "Y" ~ "sALS"),
         Case.Number = str_remove_all(Case.Number, "-"))

colnames(demographics) <- c("sample_id", "sex", "age", "clinical_diagnosis", "c9orf72_mutation", "group")

key <- key %>% 
  left_join(demographics,
            by = "sample_id")

# Prep quant data ---------------------------------------------------------

files <- list.files("myeloid_quant",
                    recursive = F,
                    full.names = T)
ids <- str_split_i(files, "/", i = 2) %>% 
  str_split_i("_", i = 1) %>%
  str_remove_all("JSB")

tabs <- list()

for (i in seq_along(files)){
  
  message("Task ", i, " out of ", length(files))
  
  hdr <- make.names(names(data.table::fread(files[i], nrows = 0)), unique = T)
  keep <- which(hdr %in% c("Object.Id", "XMin", "XMax", "YMin", "YMax") |
                  startsWith(hdr, "Region.Area") |
                  str_detect(hdr, "Average.Positive.Intensity"))
  
  tabs[[i]] <- data.table::fread(files[i], select = keep, data.table = F, col.names = hdr[keep])
  
  tabs[[i]] <- tabs[[i]] %>% 
    mutate(x = (XMax + XMin) / 2,
           y = (YMax + YMin) / 2,
           code = ids[i])
  
  tabs[[i]] <- tabs[[i]] %>% 
    mutate(Object.Id = paste0(code, "_", Object.Id)) %>% 
    dplyr::rename("object_id" = "Object.Id",
                  "area" = "Region.Area..μm..")
}

tab <- list_rbind(tabs)

colnames(tab) <- str_remove_all(colnames(tab), ".Average.Positive.Intensity") %>% 
  str_replace_all("[..]", "-") %>% 
  str_remove_all("-")

tab <- tab %>% 
  dplyr::select(-c(XMin, XMax, YMin, YMax)) %>%
  left_join(key,
            by = "code")

# Preliminary filtering on size -------------------------------------------

# tab %>% 
#   ggplot(aes(x = area)) + 
#   geom_histogram(fill = "gold",
#                  color = "black") + 
#   geom_vline(color = "midnightblue",
#              xintercept = 35) + 
#   geom_vline(color = "midnightblue",
#              xintercept = 300) +
#   theme_linedraw() + 
#   scale_x_log10() +
#   labs(x = "Square microns",
#        y = "Count") +
#   theme_linedraw() 

tab %>% 
  mutate(size_bin = case_when(area < 65 ~ "too small",
                              .default = "good")) %>% 
  group_by(size_bin) %>% 
  summarize(n = n())

size_lower <- 65
# size_upper <- 300

tab <- tab %>% 
  filter(area > size_lower)

# Prep SPE object ---------------------------------------------------------
cell_id <- tab$object_id

meta <- tab %>% 
  dplyr::select(c(object_id, area, x, y, code, sample_id, tissue, sex, age, clinical_diagnosis, c9orf72_mutation, group)) %>%
  mutate(batch = str_split_i(code, "-", i = 1))

counts <- tab %>%
  dplyr::select(-c(object_id, area, x, y, code, sample_id, tissue, sex, age, clinical_diagnosis, c9orf72_mutation, group)) %>% 
  t()

rownames(meta) <- cell_id -> colnames(counts)

sce <- SpatialExperiment(assays = list(counts = counts),
                         colData = meta,
                         spatialCoordsNames = c("x", "y"))

# Arcsinh transform -------------------------------------------------------

set.seed(1)
cof_sub <- assay(sce, "counts")[, sample(ncol(sce), min(ncol(sce), 5e4))]

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
                      exprs_values = "asinh",
                      subset_row = c("HLAA", "CD68", "CD44", "Vimentin", "CD45", "CD11c",
                                     "TMEM119", "Iba1", "Ki67", "HLADR", "iNOS", "PCNA",
                                     "GPNMB", "H2AX", "CD14", "ApoE", "CD74", "CD163", "Mac2Galectin3",
                                     "ASC", "SPP1", "C1QA", "pNRF2", "MMP9", "PSAP", "APOC1", "STING", "LYVE1"))
sce

# scater::plotPCA(sce,
#         colour_by = "group")
# 
# scater::plotPCA(sce,
#         colour_by = "batch")
# 
# scater::plotPCA(sce,
#         colour_by = "code")

# percent.var <- attr(reducedDim(sce), "percentVar")
# plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
num_pcs <- 8

# # plotReducedDim(sce, dimred = "PCA", ncomponents = 3, colour_by = "sample_id")
# 
# explain_pcs <- getExplanatoryPCs(sce,
#                                  variables = c("code",
#                                                "batch",
#                                                "group",
#                                                "sex", "age",
#                                                "area")
# )
# 
# plotExplanatoryPCs(explain_pcs/100)

# sce <- scater::runUMAP(sce,
#                dimred = "PCA",
#                ntop = num_pcs)

# scater::plotUMAP(sce,
#          colour_by = "group")
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

dittoDimPlot(sce,
             "code",
             "HARMONY_UMAP",
             split.by = "batch")

dittoDimPlot(sce, "area", "HARMONY_UMAP")

dittoDimPlot(sce, "CD74", "HARMONY_UMAP",
             assay = "asinh")

dittoPlot(sce,
          var = "Iba1",
          group.by = "code",
          plots = "vlnplot",
          assay = "counts")

dittoPlot(sce,
          var = "Iba1",
          group.by = "code",
          plots = "vlnplot",
          assay = "asinh")

dittoPlot(sce,
          var = "area",
          group.by = "code",
          plots = "vlnplot")

# DE? ---------------------------------------------------------------------

sce$group <- factor(sce$group,
                    levels = c("Control", "sALS", "C9orf72-ALS"))

de <- presto::wilcoxauc(sce, group_by = "group",
                        assay = "asinh")

de %>%
  filter(padj < 0.05)

sals <- presto::wilcoxauc(sce,
                          group_by = "group",
                          groups_use = c("Control", "sALS"),
                          assay = "asinh")

c9 <- presto::wilcoxauc(sce,
                        group_by = "group",
                        groups_use = c("Control", "C9orf72-ALS"),
                        assay = "asinh")

sals <- sals %>% 
  filter(group == "sALS") %>% 
  dplyr::select(c(feature, group, avgExpr, logFC, pval))

c9 <- c9 %>% 
  filter(group == "C9orf72-ALS") %>% 
  dplyr::select(c(feature, group, avgExpr, logFC, pval))

tab <- full_join(sals, c9,
                 by = "feature")

# padj <- p.adjust(tab$pval,
#                  method = "BH")
# 
# tab$padj <- padj
# 
# tab %>% 
#   ggplot(aes(x = padj)) + 
#   geom_histogram()
# 
# tab %>%
#   ggplot(aes(x = logFC)) + 
#   geom_histogram()

# tab %>% 
#   filter(group == "C9orf72-ALS") %>% 
#   ggplot(aes(x = logFC,
#              y = -log10(padj))) + 
#   geom_point() + 
#   theme_linedraw()

tab %>% 
  ggplot(aes(x = logFC.x,
             y = logFC.y)) + 
  geom_point() + 
  theme_linedraw()

tab %>% 
  arrange(desc(logFC.x)) %>% 
  head()

# Random plots ------------------------------------------------------------

meta %>% 
  filter(group %in% c("sALS", "Control")) %>%
  ggplot(aes(x = x,
             y = y)) + 
  geom_point(size = 0.01) + 
  facet_wrap(. ~ group,
             scales = "free",
             ncol = 1) +
  theme_void()

dittoDimPlot(sce, "CD163", assay = "asinh")
# colData(sce)
# sce
dittoPlot(sce, "TMEM119", group.by = "group",
          plots = c("vlnplot", "boxplot"),
          assay = "asinh") + 
  ylab("Normalized intensity") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
