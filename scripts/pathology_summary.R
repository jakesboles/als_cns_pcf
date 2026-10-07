library(tidyverse)
library(ggplot2)
library(janitor)
library(ggbeeswarm)
library(ggpubr)

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

# Prep Halo summary data --------------------------------------------------

dpr <- read.csv("pathology_quant/dpr_summary_data.csv")
colnames(dpr)

ptdp <- read.csv("pathology_quant/ptdp_summary_data.csv")

ptdp <- ptdp[, c("Image.Tag", colnames(ptdp)[startsWith(colnames(ptdp), "ptdp")])] %>% 
  na.omit()
dpr <- dpr[, c("Image.Tag", colnames(dpr)[startsWith(colnames(dpr), "dpr")])] %>% 
  na.omit()

head(ptdp)
head(dpr)

colnames(ptdp) <- c("image", "analyzed_area", "bg_area", "feature_area") -> colnames(dpr)
ptdp$feature <- "ptdp"
dpr$feature <- "dpr"

df <- rbind(ptdp, dpr)

head(df)

df <- df %>% 
  mutate(code = str_split_i(image, "_", i = 1) %>% 
           str_remove_all("JSB"),
         feature_pct = (feature_area / analyzed_area) * 100)

# add demographics

df <- df %>%
  left_join(key,
            by = "code")

{roi <- "ptdp"
  
  df %>% 
    filter(feature == roi) %>%
    mutate(group = factor(group,
                          levels = c("Control", "sALS", "C9orf72-ALS")),
           tissue = factor(tissue,
                           levels = c("mcx", "sc"),
                           labels = c("Motor cortex", "Spinal cord"))) %>%
    ggplot(aes(x = group,
               y = feature_pct,
               fill = group)) +
    geom_quasirandom(shape = 21,
                     show.legend = F,
                     size = 5,
                     alpha = 0.7) +
    stat_summary(fun = mean,
                 geom = "crossbar",
                 width = 0.7,
                 show.legend = F) + 
    stat_summary(fun.data = mean_se,
                 geom = "errorbar",
                 width = 0.5, 
                 linewidth = 1.2,
                 show.legend = F) +
    scale_fill_manual(values = c("#b8b0a8", "#CC00FF", "#0CAA00")) +
    geom_pwc(method = "wilcox_test",
             hide.ns = F,
             # step.increase = step_increase,
             label = "{p.adj.signif}") +
    # ggtitle(roi) +
    # ggtitle("DPR inclusions") +
    ggtitle("pTDP-43 inclusions") +
    labs(y = "% area") +
  facet_wrap(. ~ tissue) +
  theme_linedraw(base_size = 16) + 
    theme(axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  }
