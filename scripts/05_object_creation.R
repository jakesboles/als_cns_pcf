dyn.load("/hpc/software/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libxml2-2.9.10-sbmpowvpmrw3sytyglnc3y6stcfbkg4o/lib/libxml2.so.2", local = F)
suppressPackageStartupMessages({
  library(tidyverse)
  library(SpatialExperiment)
})

setwd("/projects/b1169/boles/als_cns_pcf")

data_dir <- "data/05_object_creation/"
results_dir <- "results/05_object_creation/"
dir.create(data_dir, showWarnings = F, recursive = T)
dir.create(results_dir, showWarnings = F, recursive = T)

# Settings ----------------------------------------------------------------

tissues_use <- "mcx"

channel_rename <- c("Mac2/Galectin-3" = "Galectin-3")

# channels to leave out of the object, by name after renaming, e.g. c("Ki67", "H2A.X")
channels_drop <- c("DAPI", "Tau (phospho T181)", "Tau (phospho S396)",
                           "EPR2488", "Tau (phospho T205)", "A-beta")

canonical_markers <- list(
  astrocyte = c("Vimentin", "GFAP"),
  cd4_tcell = c("CD4", "CD3e"),
  cd8_tcell = c("CD8", "CD3e"),
  cortical_neuron = c("NeuN"),
  myeloid = c("Iba-1"),
  oligodendrocyte = c("Olig-2")
)

# checked reasonable size limits in Qupath 11/5/26
size_limits <- tribble(
  ~cell_type,         ~min_area, ~max_area,
  "astrocyte",       40,        440,
  "cd4_tcell",        10,        Inf,
  "cd8_tcell",        10,        Inf,
  "cortical_neuron", 40,        Inf,
  "myeloid",          38,        Inf,
  "oligodendrocyte",  8,        200
)

cells_per_section <- 5000
cutoff_k <- 2

# extra cell types expected to express a canonical marker, left out of its background reference,
# e.g. list(CD4 = "myeloid", Vimentin = "myeloid")
marker_extra_positive <- list()

# manual cofactor overrides, e.g. c(`Iba-1` = 120); thresholds are tuned below
cofactor_manual <- c()

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

for (i in seq_along(channels)){
  if (anyDuplicated(channels[[i]]) > 0) { 
    message(paste0("Duplicated channel in ", names(channels)[i]))
  }
}

features <- channels[[1]]
mismatch <- purrr::keep(channels, ~ !setequal(.x, features))
if (length(mismatch)) {
  walk(names(mismatch), function(cd) {
    message(cd, ": missing ", paste(setdiff(features, mismatch[[cd]]), collapse = ", "),
            " | extra ", paste(setdiff(mismatch[[cd]], features), collapse = ", "))
  })
  stop("Channel names differ between sections; extend channel_rename")
}

unknown_drop <- setdiff(channels_drop, features)
if (length(unknown_drop)) stop("channels_drop not found in channels: ", paste(unknown_drop, collapse = ", "))
features <- setdiff(features, channels_drop)
message(length(features), " channels kept, ", length(channels_drop), " dropped")

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
  if (anyDuplicated(names(d1)) || !all(c(features, channels_drop) %in% names(d1)))
    stop(cd, ": 01 parquet columns don't match its _channels.csv (duplicates or missing channels); rerun 01")
  d1 <- d1 %>% dplyr::select(-any_of(channels_drop))
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
cof_types <- spe$cell_type[cof_idx]

canon <- unique(unlist(canonical_markers))
bad_extra <- setdiff(names(marker_extra_positive), canon)
if (length(bad_extra)) stop("marker_extra_positive names aren't canonical markers: ", paste(bad_extra, collapse = ", "))
bad_types <- setdiff(unlist(marker_extra_positive), names(canonical_markers))
if (length(bad_types)) stop("marker_extra_positive cell types unknown: ", paste(bad_types, collapse = ", "))

pos_types <- map(set_names(canon), ~ c(names(purrr::keep(canonical_markers, function(m) .x %in% m)),
                                        marker_extra_positive[[.x]]))
ref_cells <- map(pos_types, ~ !(cof_types %in% .x))

bg_peak <- function(v) {
  d <- density(log(v[!is.na(v) & v > 0]), n = 512)
  pk <- which(diff(sign(diff(d$y))) == -2) + 1
  pk <- pk[d$y[pk] > 0.05 * max(d$y)]
  if (!length(pk)) pk <- which.max(d$y)
  exp(d$x[min(pk)])
}

cofactor <- apply(cof_sub, 1, bg_peak)
for (mk in canon) cofactor[mk] <- bg_peak(cof_sub[mk, ref_cells[[mk]]])
if (length(cofactor_manual)) cofactor[names(cofactor_manual)] <- cofactor_manual
stopifnot(identical(names(cofactor), rownames(spe)))

# Automatic background cutoffs (asinh scale) ------------------------------

a0 <- asinh(1)
asinh_sub <- asinh(cof_sub / cofactor)
left_sd <- function(a) {
  a <- a[!is.na(a) & a < a0]
  if (length(a)) sqrt(mean((a - a0)^2)) else 0
}
bg_sd <- apply(asinh_sub, 1, left_sd)
for (mk in canon) bg_sd[mk] <- left_sd(asinh_sub[mk, ref_cells[[mk]]])
threshold_auto <- cofactor * sinh(a0 + cutoff_k * bg_sd)

p <- as_tibble(t(cof_sub)) %>%
  pivot_longer(everything(), names_to = "marker") %>%
  filter(value > 0) %>%
  ggplot(aes(x = value)) +
  geom_density() +
  geom_vline(data = enframe(cofactor, "marker", "cofactor"), aes(xintercept = cofactor), color = "firebrick") +
  scale_x_log10() +
  facet_wrap(. ~ marker, scales = "free") +
  theme_linedraw(base_size = 8)
ggsave(paste0(results_dir, "cofactor_density.png"), p, width = 20, height = 16, dpi = 200)

rm(cof_sub, asinh_sub)
gc()

# Tune canonical marker thresholds ----------------------------------------

set.seed(1)
plot_idx <- split(seq_len(ncol(spe)), spe$cell_type) %>%
  map(~ if (length(.x) > 2e4) sample(.x, 2e4) else .x) %>%
  unlist(use.names = F)

lineage_long <- map(canon, function(mk) {
  tibble(marker = mk,
         value = assay(spe, "counts")[mk, plot_idx],
         cell_type = spe$cell_type[plot_idx],
         expected = spe$cell_type[plot_idx] %in% pos_types[[mk]])
}) %>%
  list_rbind() %>%
  filter(value > 0)

plot_thresholds <- function(thr) {
  lines <- tibble(marker = canon, cofactor = cofactor[canon], auto = threshold_auto[canon], threshold = thr[canon])
  lineage_long %>%
    ggplot(aes(x = value, fill = expected)) +
    geom_density(alpha = 0.6) +
    geom_vline(data = lines, aes(xintercept = cofactor), color = "firebrick") +
    geom_vline(data = lines, aes(xintercept = auto), color = "gray40", linetype = 3) +
    geom_vline(data = lines, aes(xintercept = threshold), color = "dodgerblue", linetype = 2) +
    scale_x_log10() +
    scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "gray70"),
                      labels = c(`TRUE` = "expected positive", `FALSE` = "background reference")) +
    facet_grid(cell_type ~ marker, scales = "free") +
    labs(x = "raw intensity", fill = NULL) +
    theme_linedraw(base_size = 8)
}

lineage_negative <- function(thr) {
  neg <- rep(FALSE, ncol(spe))
  for (ct in names(canonical_markers)) {
    sel <- which(spe$cell_type == ct)
    mk <- canonical_markers[[ct]]
    neg[sel] <- colSums(assay(spe, "counts")[mk, sel, drop = F] >= thr[mk]) == 0
  }
  neg
}

preview_flags <- function(thr) {
  tibble(cell_type = spe$cell_type, group = spe$group, flagged = lineage_negative(thr)) %>%
    group_by(cell_type, group) %>%
    summarize(pct_flagged = round(100 * mean(flagged), 1), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = pct_flagged)
}

# raw intensity read off the plot's x-axis, e.g. c(`Iba-1` = 800, NeuN = 150)
# edit and re-run from here through preview_flags() until the blue lines sit where you want
threshold_manual <- c('CD3e' = (cofactor[['CD3e']] * threshold_auto[['CD3e']])^(1/2),
                            'CD4' = threshold_auto[['CD4']],
                            'CD8' = threshold_auto[['CD8']],
                            'GFAP' = cofactor[["GFAP"]],
                            'Iba-1' = (cofactor[['Iba-1']] * threshold_auto[["Iba-1"]])^(1/2),
                            'NeuN' = cofactor[['NeuN']],
                            'Olig-2' = threshold_auto[["Olig-2"]],
                            'Vimentin' = cofactor[["Vimentin"]])

bad_thr <- setdiff(names(threshold_manual), canon)
if (length(bad_thr)) stop("threshold_manual names aren't canonical markers: ", paste(bad_thr, collapse = ", "))
threshold <- threshold_auto
if (length(threshold_manual)) threshold[names(threshold_manual)] <- threshold_manual

plot_thresholds(threshold)
preview_flags(threshold)

# Save thresholds and transform -------------------------------------------

ggsave(paste0(results_dir, "canonical_marker_by_cell_type.png"), plot_thresholds(threshold),
       width = 16, height = 10, dpi = 200)

marker_table <- tibble(marker = names(cofactor), cofactor = cofactor, bg_sd = bg_sd,
                       threshold_auto = threshold_auto, threshold = threshold,
                       threshold_asinh = asinh(threshold / cofactor),
                       manual = names(cofactor) %in% names(threshold_manual),
                       background_from = if_else(names(cofactor) %in% canon,
                                                 "other cell types", "all cells"))
write.csv(marker_table, paste0(results_dir, "marker_cofactors_thresholds.csv"), row.names = F)

assay(spe, "asinh") <- asinh(assay(spe, "counts") / cofactor)

# Z normalization (Hickey 2021; SPACEc 2025): per marker over all cells and sections, raw intensities, no cap
z_mean <- rowMeans(assay(spe, "counts"), na.rm = T)
z_sd <- matrixStats::rowSds(assay(spe, "counts"), na.rm = T)
names(z_sd) <- rownames(spe)
stopifnot(identical(names(z_mean), rownames(spe)), all(z_sd > 0))
assay(spe, "z") <- (assay(spe, "counts") - z_mean) / z_sd
write.csv(tibble(marker = rownames(spe), z_mean = z_mean, z_sd = z_sd),
          paste0(results_dir, "marker_z_params.csv"), row.names = F)

# Flags -------------------------------------------------------------------

limits <- size_limits[match(spe$cell_type, size_limits$cell_type), ]
spe$flag_size <- spe$area_um2 <= limits$min_area | spe$area_um2 >= limits$max_area

spe$flag_region <- is.na(spe$region) | spe$region %in% c("unassigned", "multiple")

spe$flag_lineage <- lineage_negative(threshold)

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

table(spe$cell_type, spe$flag_reason)

# Save --------------------------------------------------------------------

saveRDS(spe, paste0(data_dir, "spe.rds"))
spe

# Export flags for 06 -----------------------------------------------------

arrow::write_parquet(as.data.frame(colData(spe))[, c("code", "cell_id", "cell_type", "area_um2",
                                                     flag_cols, "keep", "flag_reason")],
                     paste0(data_dir, "cell_flags.parquet"))
