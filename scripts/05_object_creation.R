dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
library(tidyverse)
library(SpatialExperiment)

setwd("/projects/b1169/boles/als_cns_pcf")

data_dir <- "data/05_object_creation/"
results_dir <- "results/05_object_creation/"
dir.create(data_dir, showWarnings = F, recursive = T)
dir.create(results_dir, showWarnings = F, recursive = T)

# Settings ----------------------------------------------------------------

tissues_use <- "mcx"

channel_rename <- c("Mac2/Galectin-3" = "Galectin-3")

canonical_markers <- list(
  astrocytes = c("Vimentin", "GFAP"),
  cd4_tcell = c("CD4", "CD3e"),
  cd8_tcell = c("CD8", "CD3e"),
  cortical_neurons = c("NeuN"),
  myeloid = c("Iba-1"),
  oligodendrocyte = c("Olig-2")
)

size_limits <- tribble(
  ~cell_type,         ~min_area, ~max_area,
  "astrocytes",       40,        250,
  "cd4_tcell",        10,        Inf,
  "cd8_tcell",        10,        Inf,
  "cortical_neurons", 40,        400,
  "myeloid",          40,        400,
  "oligodendrocyte",  10,        60
)

cells_per_section <- 5000
cutoff_k <- 3

# manual overrides, e.g. c(`Iba-1` = 120)
cofactor_manual <- c()
cutoff_manual <- c()

# Sample key and demographics ---------------------------------------------

key <- read.csv("sample_key.csv", fileEncoding = "UTF-8-BOM") %>%
  dplyr::select(code, sample_id, tissue) %>%
  mutate(across(everything(), str_trim))

demographics <- read.csv("target_als_demographics_compiled.csv") %>%
  dplyr::select(Case.Number, Sex, Age.at.Death, Clinical.Diagnosis, C9orf72.mutation) %>%
  mutate(group = case_when(C9orf72.mutation == "Y" ~ "C9orf72-ALS",
                           Clinical.Diagnosis == "Control" ~ "Control",
                           .default = "sALS"),
         Case.Number = str_remove_all(Case.Number, "-")) %>%
  dplyr::rename(sample_id = Case.Number, sex = Sex, age = Age.at.Death,
                clinical_diagnosis = Clinical.Diagnosis, c9orf72_mutation = C9orf72.mutation)

key <- key %>%
  left_join(demographics, by = "sample_id") %>%
  mutate(batch = str_split_i(code, "-", i = 1),
         group = factor(group, levels = c("Control", "sALS", "C9orf72-ALS")))

codes <- list.files("data/01_quantify_cells", pattern = "\\.parquet$") %>%
  str_remove("\\.parquet$")
codes <- intersect(codes, key$code[key$tissue %in% tissues_use])
message(length(codes), " sections: ", paste(codes, collapse = ", "))

# Channel names -----------------------------------------------------------

rename_channels <- function(x) {
  i <- x %in% names(channel_rename)
  x[i] <- channel_rename[x[i]]
  unname(x)
}

channels <- map(set_names(codes), function(cd) {
  rename_channels(read.csv(paste0("results/01_quantify_cells/", cd, "_channels.csv"))$channel)
})

features <- channels[[1]]
mismatch <- purrr::keep(channels, ~ !setequal(.x, features))
if (length(mismatch)) {
  walk(names(mismatch), function(cd) {
    message(cd, ": missing ", paste(setdiff(features, mismatch[[cd]]), collapse = ", "),
            " | extra ", paste(setdiff(mismatch[[cd]], features), collapse = ", "))
  })
  stop("Channel names differ between sections; extend channel_rename")
}

missing_markers <- setdiff(unlist(canonical_markers), features)
if (length(missing_markers)) stop("Canonical markers not in channels: ", paste(missing_markers, collapse = ", "))

# Load 01-04 --------------------------------------------------------------

read_step <- function(step, cd) {
  f <- paste0("data/", step, "/", cd, ".parquet")
  if (!file.exists(f)) stop("Missing ", f, "; run ", step, " for ", cd)
  arrow::read_parquet(f)
}

pq <- map(codes, function(cd) {
  message("Loading ", cd)
  d1 <- read_step("01_quantify_cells", cd)
  names(d1) <- rename_channels(names(d1))
  d2 <- read_step("02_annotate_regions", cd) %>% dplyr::select(cell_id, in_gm, in_wm, region)
  d3 <- read_step("03_annotate_pathology", cd) %>% dplyr::select(cell_id, any_of(c("has_dpr", "has_ptdp")))
  d4 <- read_step("04_feature_distances", cd) %>% dplyr::select(-code)
  out <- d1 %>%
    left_join(d2, by = "cell_id") %>%
    left_join(d3, by = "cell_id") %>%
    left_join(d4, by = "cell_id")
  stopifnot(nrow(out) == nrow(d1))
  out
}) %>%
  list_rbind()

pq <- pq %>%
  left_join(key, by = "code")

unknown_types <- setdiff(unique(pq$cell_type), size_limits$cell_type)
if (length(unknown_types)) stop("No size limits / markers for cell types: ", paste(unknown_types, collapse = ", "))

# Build SPE ---------------------------------------------------------------

meta <- pq %>%
  dplyr::select(-all_of(features)) %>%
  as.data.frame()
rownames(meta) <- meta$cell_id

counts <- t(as.matrix(pq[, features]))
colnames(counts) <- meta$cell_id

rm(pq)
gc()

spe <- SpatialExperiment(assays = list(counts = counts),
                         colData = meta,
                         spatialCoordsNames = c("x", "y"))
rm(counts, meta)
gc()

# Cofactors (background peak) ---------------------------------------------

set.seed(1)
cof_idx <- split(seq_len(ncol(spe)), spe$code) %>%
  map(~ if (length(.x) > cells_per_section) sample(.x, cells_per_section) else .x) %>%
  unlist(use.names = F)
cof_sub <- assay(spe, "counts")[, cof_idx]

bg_peak <- function(v) {
  d <- density(log(v[!is.na(v) & v > 0]), n = 512)
  pk <- which(diff(sign(diff(d$y))) == -2) + 1
  pk <- pk[d$y[pk] > 0.05 * max(d$y)]
  if (!length(pk)) pk <- which.max(d$y)
  exp(d$x[min(pk)])
}

cofactor <- apply(cof_sub, 1, bg_peak)
if (length(cofactor_manual)) cofactor[names(cofactor_manual)] <- cofactor_manual
stopifnot(identical(names(cofactor), rownames(spe)))

assay(spe, "asinh") <- asinh(assay(spe, "counts") / cofactor)

# Background cutoffs (asinh scale) ----------------------------------------

a0 <- asinh(1)
asinh_sub <- asinh(cof_sub / cofactor)
bg_sd <- apply(asinh_sub, 1, function(a) {
  a <- a[!is.na(a) & a < a0]
  if (length(a)) sqrt(mean((a - a0)^2)) else 0
})
cutoff <- a0 + cutoff_k * bg_sd
if (length(cutoff_manual)) cutoff[names(cutoff_manual)] <- cutoff_manual

marker_table <- tibble(marker = names(cofactor), cofactor = cofactor, bg_sd = bg_sd,
                       cutoff_asinh = cutoff, cutoff_raw = cofactor * sinh(cutoff))
write.csv(marker_table, paste0(results_dir, "marker_cofactors_cutoffs.csv"), row.names = F)

p <- as_tibble(t(cof_sub)) %>%
  pivot_longer(everything(), names_to = "marker") %>%
  filter(value > 0) %>%
  ggplot(aes(x = value)) +
  geom_density() +
  geom_vline(data = marker_table, aes(xintercept = cofactor), color = "firebrick") +
  geom_vline(data = marker_table, aes(xintercept = cutoff_raw), color = "dodgerblue", linetype = 2) +
  scale_x_log10() +
  facet_wrap(. ~ marker, scales = "free") +
  theme_linedraw(base_size = 8)
ggsave(paste0(results_dir, "cofactor_density.png"), p, width = 20, height = 16, dpi = 200)

rm(cof_sub, asinh_sub)
gc()

# Flags -------------------------------------------------------------------

limits <- size_limits[match(spe$cell_type, size_limits$cell_type), ]
spe$flag_size <- spe$area_um2 <= limits$min_area | spe$area_um2 >= limits$max_area

spe$flag_region <- is.na(spe$region) | spe$region %in% c("unassigned", "multiple")

flag_lineage <- rep(FALSE, ncol(spe))
for (ct in names(canonical_markers)) {
  sel <- which(spe$cell_type == ct)
  mk <- canonical_markers[[ct]]
  pos <- assay(spe, "asinh")[mk, sel, drop = F] >= cutoff[mk]
  flag_lineage[sel] <- colSums(pos) == 0
}
spe$flag_lineage <- flag_lineage

flag_cols <- c("flag_size", "flag_region", "flag_lineage")
flags <- as.matrix(as.data.frame(colData(spe)[, flag_cols]))
spe$keep <- rowSums(flags) == 0
reason <- rep("", ncol(spe))
for (f in flag_cols) {
  hit <- flags[, f] %in% TRUE
  lab <- str_remove(f, "flag_")
  reason[hit] <- if_else(reason[hit] == "", lab, paste(reason[hit], lab, sep = ";"))
}
spe$flag_reason <- reason
rm(flags, reason)

# Flag QC -----------------------------------------------------------------

flag_summary <- as.data.frame(colData(spe)) %>%
  group_by(code, sample_id, group, cell_type) %>%
  summarize(n = n(),
            across(all_of(flag_cols), ~ mean(.x)),
            kept = mean(keep),
            .groups = "drop")
write.csv(flag_summary, paste0(results_dir, "flag_summary.csv"), row.names = F)

p <- flag_summary %>%
  pivot_longer(all_of(flag_cols), names_to = "flag", values_to = "frac") %>%
  ggplot(aes(x = group, y = frac, color = group)) +
  geom_jitter(width = 0.15, height = 0) +
  facet_grid(flag ~ cell_type, scales = "free_y") +
  labs(y = "Fraction of cells flagged (per section)") +
  theme_linedraw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(paste0(results_dir, "flag_fraction_by_group.png"), p, width = 12, height = 9, dpi = 200)

set.seed(1)
plot_idx <- split(seq_len(ncol(spe)), spe$cell_type) %>%
  map(~ if (length(.x) > 2e4) sample(.x, 2e4) else .x) %>%
  unlist(use.names = F)

lineage_long <- map(unique(unlist(canonical_markers)), function(mk) {
  own_types <- names(purrr::keep(canonical_markers, ~ mk %in% .x))
  tibble(marker = mk,
         value = assay(spe, "asinh")[mk, plot_idx],
         own = spe$cell_type[plot_idx] %in% own_types)
}) %>%
  list_rbind()

p <- lineage_long %>%
  ggplot(aes(x = value, color = own)) +
  geom_density() +
  geom_vline(data = filter(marker_table, marker %in% unlist(canonical_markers)),
             aes(xintercept = cutoff_asinh), linetype = 2) +
  facet_wrap(. ~ marker, scales = "free") +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "gray50"),
                     labels = c(`TRUE` = "own cell type", `FALSE` = "other cell types")) +
  labs(x = "asinh intensity", color = NULL) +
  theme_linedraw(base_size = 9)
ggsave(paste0(results_dir, "canonical_marker_cutoffs.png"), p, width = 10, height = 6, dpi = 200)

table(spe$cell_type, spe$flag_reason)

# Save --------------------------------------------------------------------

saveRDS(spe, paste0(data_dir, "spe.rds"))
spe
