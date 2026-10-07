library(tidyverse)
library(ggplot2)
library(janitor)
library(ggbeeswarm)
library(ggpubr)
library(SingleCellExperiment)
library(scuttle)
library(scater)
library(presto)
library(dittoSeq)

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

files <- list.files("pathology_quant",
                    pattern = "ptdp",
                    recursive = F,
                    full.names = T)
ids <- str_split_i(files, "/", i = 2) %>% 
  str_split_i("_", i = 1) %>%
  str_remove_all("JSB")

tabs <- list()

for (i in seq_along(files)){
  
  message("Task ", i, " out of ", length(files))
  
  tabs[[i]] <- read.csv(files[i])
  
  tabs[[i]] <- tabs[[i]] %>% 
    mutate(x = (XMax + XMin) / 2,
           y = (YMax + YMin) / 2,
           code = ids[i])
  
  feature_cols <- colnames(tabs[[i]])[str_detect(colnames(tabs[[i]]), "Average.Positive.Intensity")]
  other_cols <- c("Object.Id", "Region.Area..μm..", "x", "y", "code")
  
  tabs[[i]] <- tabs[[i]][, c(other_cols, feature_cols)] %>%
    mutate(Object.Id = paste0(code, "_", Object.Id)) %>% 
    dplyr::rename("object_id" = "Object.Id",
                  "area" = "Region.Area..μm..")
}

tab <- list_rbind(tabs)

feature_col_names <- str_remove_all(feature_cols, ".Average.Positive.Intensity") %>% 
  str_replace_all("[..]", "-") %>%
  str_remove_all("-")

colnames(tab)

colnames(tab) <- c("object.id", "area", "x", "y", "code", feature_col_names)

tab <- tab %>% 
  left_join(key,
            by = "code")

# Preliminary filtering on size -------------------------------------------

# tab %>% 
#   ggplot(aes(x = area)) + 
#   geom_histogram(binwidth = 10) + 
#   theme_linedraw() + 
#   scale_x_continuous(expand = c(0, 0))
# 
# tab %>% 
#   mutate(size_bin = case_when(area < 25 ~ "too small",
#                               area >= 20 & area < 200 ~ "good",
#                               area >= 200 ~ "too big")) %>% 
#   group_by(size_bin) %>% 
#   summarize(n = n())
# 
# size_lower <- 25
# size_upper <- 200
# 
# tab <- tab %>% 
#   filter(area > size_lower & area < size_upper)

# Prep SCE object ---------------------------------------------------------
cell_id <- tab$object.id

meta <- tab %>% 
  dplyr::select(c(object.id, area, x, y, code, sample_id, tissue, sex, age, clinical_diagnosis, c9orf72_mutation, group))

counts <- tab %>%
  dplyr::select(-c(object.id, area, x, y, code, sample_id, tissue, sex, age, clinical_diagnosis, c9orf72_mutation, group)) %>% 
  t()

rownames(meta) <- cell_id -> colnames(counts)

sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = meta)

sce <- scuttle::logNormCounts(sce)
sce

sce <- scater::runPCA(sce)
sce

plotPCA(sce,
        colour_by = "group")

percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
num_pcs <- 6

plotReducedDim(sce, dimred = "PCA", ncomponents = 3, colour_by = "sample_id")
# 
# explain_pcs <- getExplanatoryPCs(sce,
#                                  variables = c("code")
# )
# 
# plotExplanatoryPCs(explain_pcs/100)

# Extract PCA cell embeddings
pca_coords <- reducedDim(sce, "PCA")

# Extract the feature loadings matrix attached as an attribute
pca_loadings <- attr(pca_coords, "rotation") 
head(pca_loadings)

pca_loadings %>% 
  as.data.frame() %>%
  rownames_to_column("feature") %>% 
  ggplot() + 
  geom_segment(aes(xend = PC1,
                   yend = PC2,
                   x = 0, y = 0))

sce <- runUMAP(sce)

plotUMAP(sce, colour_by = "code")

# DE? ---------------------------------------------------------------------

de <- presto::wilcoxauc(sce, group_by = "group")

de %>%
  filter(abs(logFC) > 0.5)

sals <- de %>% 
  filter(group == "sALS") %>%
  dplyr::select(c(feature, group, logFC))

c9 <- de %>%
  filter(group == "C9orf72-ALS") %>%
  dplyr::select(c(feature, group, logFC))

tab <- full_join(sals,  c9,
                 by = "feature")

tab %>% 
  ggplot(aes(x = logFC.x,
             y = logFC.y)) + 
  geom_point()

# Random plots ------------------------------------------------------------

sce$group <- factor(sce$group,
                    levels = c("Control", "sALS", "C9orf72-ALS"))

meta %>% 
  filter(group %in% c("sALS", "Control")) %>%
  ggplot(aes(x = x,
             y = y)) + 
  geom_point(size = 0.01) + 
  facet_wrap(. ~ group,
             scales = "free",
             ncol = 1) +
  theme_void()

dittoDimPlot(sce, "CD163")
# colData(sce)
# sce
dittoPlot(sce, "TMEM119", group.by = "group",
          plots = c("vlnplot", "boxplot")) + 
  ylab("Normalized intensity") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
