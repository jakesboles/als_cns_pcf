suppressPackageStartupMessages({ library(tidyverse); library(lme4) })
setwd("/projects/b1169/boles/als_cns_pcf")

source("/projects/b1169/boles/als_cns_scrnaseq/r_scripts/plot_theme.R")

out <- "results/13_inclusion_neighborhood/"; dir.create(out, showWarnings = F, recursive = T)

rings  <- c("0", "0-5", "5-15", "15-30", "30-60", "60-90", "90-120", "120-150", ">150")
b <- read_csv("data/13_inclusion_neighborhood/band_counts.csv", show_col_types = F) %>%
  filter(group %in% c("sALS", "C9orf72-ALS")) %>%
  mutate(ring = factor(band, rings), group = factor(group, c("sALS", "C9orf72-ALS")),
         label = recode(label, ptdp = "pTDP-43", dpr = "DPR"),
         log2_oe = log2((n_cells + 0.5) / (exp_label + 0.5)))

# 1. QC
inc <- arrow::read_parquet("data/13_inclusion_neighborhood/inclusion_ring_counts.parquet")
inc %>% distinct(code, group, label, inclusion) %>% count(group, code, label) %>% print(n = 50)
b %>% distinct(code, label, region, ring, n_band) %>% filter(n_band < 20) %>% print(n = 50)

# 2. profiles (rings within 100 um; >100 is the reference)
prof <- b %>% filter(ring != ">150")
p <- ggplot(prof, aes(ring, log2_oe, colour = group)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(group = code), alpha = 0.3) +
  stat_summary(aes(group = group), fun = mean, geom = "line", linewidth = 1.2) +
  facet_grid(label + region ~ cell_type2) +
  labs(x = "distance to nearest inclusion (µm, centroid to edge)", y = "log2 observed / expected") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(out, "enrichment_profiles.png"), p, width = 14, height = 8, dpi = 300)

prof %>% 
  filter(cell_type2 == "Microglia" & region == "GM" & 
           ring != "0" & label == "pTDP-43") %>%
  ggplot(aes(ring, log2_oe, colour = group)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(group = code), alpha = 0.3) +
  stat_summary(aes(group = group), fun = mean, geom = "line", linewidth = 1.2) +
  # facet_grid(label + region ~ cell_type2) +
  labs(x = "µm, centroid to edge", y = "log2(observed / expected)") +
  scale_color_manual(values = c("#0CAA00", "#CC00FF")) +
  ggtitle("Microglial distance to nearest\nTDP-43 inclusion") +
  theme_bw() + 
  theme_als_project(base_size = 9) +
  rotate_x_text(30) +
  theme(legend.title = element_blank(),
        legend.box.spacing = unit(0.1, "cm"),
        legend.position = "bottom")
ggsave(filename = "results/sfn_figs/microglia_ptdp_distance.png",
       units = "in", dpi = 600,
       height = 4.5, width = 3.25)

prof %>% 
  filter(cell_type2 == "Microglia" & region == "GM" & ring != "0" & 
           group == "C9orf72-ALS") %>%
  mutate(code2 = paste0(code, "_", label)) %>%
  ggplot(aes(ring, log2_oe, colour = label)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(group = code2), alpha = 0.3) +
  stat_summary(aes(group = label), fun = mean, geom = "line", linewidth = 1.2) +
  # facet_grid(label + region ~ cell_type2) +
  labs(x = "µm, centroid to edge", y = "log2(observed / expected)") +
  # scale_color_manual(values = c("#0CAA00", "#CC00FF")) +
  ggtitle("Microglial distance to nearest\ninclusion in C9orf72-ALS") +
  scale_color_manual(values = c("red", "blue")) +
  theme_bw() + 
  theme_als_project(base_size = 9) +
  rotate_x_text(30) +
  theme(legend.title = element_blank(),
        legend.box.spacing = unit(0.1, "cm"),
        legend.position = "bottom")
ggsave(filename = "results/sfn_figs/microglia_c9_path_distance.png",
       units = "in", dpi = 600,
       height = 4.5, width = 3.25)

# 3. sALS vs C9 per label x region x cell type: Poisson GLMM on counts, offset = random-relabelling expectation
fit_one <- function(d, term) {
  d <- filter(d, exp_label > 0)
  m <- try(glmer(as.formula(paste("n_cells ~", term, "* ring + offset(log(exp_label)) + (1 | code)")),
                 family = poisson, data = d, control = glmerControl(optimizer = "bobyqa")), silent = T)
  if (inherits(m, "try-error")) return(tibble())
  m0 <- update(m, as.formula(paste(". ~ . -", paste0(term, ":ring"))))
  co <- summary(m)$coefficients
  tibble(main_effect = co[grep(paste0("^", term, "[^:]*$"), rownames(co)), "Estimate"][1] / log(2),
         p_profile = anova(m0, m)$`Pr(>Chisq)`[2])
}
tests <- prof %>% group_by(label, region, cell_type2) %>%
  group_modify(~ fit_one(.x, "group")) %>% ungroup() %>%
  mutate(padj_profile = p.adjust(p_profile, "BH"))
write_csv(tests, paste0(out, "group_tests.csv"))   # main_effect = log2 C9/sALS at ring "0"; p_profile = group:ring

# 4. DPR vs pTDP within C9 sections (paired within section)
c9 <- prof %>% filter(group == "C9orf72-ALS") %>% mutate(label = factor(label, c("pTDP-43", "DPR")))
p <- c9 %>% group_by(code, label, region, ring, cell_type2) %>% summarise(log2_oe = mean(log2_oe), .groups = "drop") %>%
  ggplot(aes(label, log2_oe, group = code)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_line(alpha = 0.5) + geom_point(size = 0.8) + facet_grid(region + cell_type2 ~ ring) +
  theme_bw(base_size = 8)
ggsave(paste0(out, "dpr_vs_ptdp_c9.png"), p, width = 10, height = 14, dpi = 300)
label_tests <- c9 %>% group_by(region, cell_type2) %>% group_modify(~ fit_one(.x, "label")) %>% ungroup() %>%
  mutate(padj_profile = p.adjust(p_profile, "BH"))
write_csv(label_tests, paste0(out, "dpr_vs_ptdp_tests.csv"))

c9 %>% 
  filter(cell_type2 == "Microglia" & region == "GM" & ring != "0") %>%
  group_by(code, label, region, ring, cell_type2) %>% 
  summarise(log2_oe = mean(log2_oe), .groups = "drop") %>%
  ggplot(aes(label, log2_oe, group = code)) + 
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(alpha = 0.5) + 
  geom_point(size = 0.8) + 
  facet_grid(region + cell_type2 ~ ring) +
  theme_bw(base_size = 8)

# 5. summary heatmap: C9 - sALS mean log2 o/e
hm <- prof %>% group_by(label, region, cell_type2, ring, group) %>% summarise(m = mean(log2_oe), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = m) %>% mutate(diff = `C9orf72-ALS` - sALS)
p <- ggplot(hm, aes(ring, cell_type2, fill = diff)) + geom_tile() +
  scale_fill_gradient2(low = "dodgerblue", high = "firebrick", name = "C9 − sALS\nlog2 o/e") +
  facet_grid(region ~ label) + theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(out, "group_difference_heatmap.png"), p, width = 8, height = 6, dpi = 300)

library(glmmTMB)
fit <- glmmTMB(n_cells ~ group * band + offset(log(exp_label)) + (1 | code),
               family = nbinom2, data = b)          # or poisson with + (1 | code:band)

# simplest honest test: one value per section, inner rings pooled
b |> filter(label == "pTDP-43", band %in% c("0", "0-5")) |>
  group_by(code, group, cell_type2) |>
  summarise(o = sum(n_cells), e = sum(exp_label), .groups = "drop") |>
  mutate(l2oe = log2((o + 0.5) / (e + 0.5))) |>
  group_by(cell_type2) |>
  summarise(diff = mean(l2oe[group == "C9orf72-ALS"]) - mean(l2oe[group == "sALS"]),
            p = wilcox.test(l2oe ~ group)$p.value)

b %>% 
  filter(cell_type2 == "Microglia" & 
           band == "0")

dittoScatterPlot(spe,
                 x.var = "GPNMB",
                 y.var = "Galectin-3")
