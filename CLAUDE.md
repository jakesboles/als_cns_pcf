# als_cns_pcf

Working notes for Claude Code sessions on this repo: observations about the data, troubleshooting history, rules of engagement, and current understanding of project intent. Running memory, not a tutorial. Update it as the project moves; don't let it go stale.

## Rules of engagement

1. **Keep this file current.** Observations about the data, troubleshooting steps worth remembering, interaction rules, and project intent all live here.
2. **Never edit code directly.** Every change (this file included) goes through a PR on the session's designated feature branch for the user to accept or reject. No need to monitor PRs; the user reports their status (may change later if GitHub-based communication becomes easier).
3. **Minimal code comments.** Short steering comments only. Explanations, rationale, and caveats go here, not in scripts. (Note: this is the opposite of the sister Visium repo's heavily-commented style; don't port its comment density.)
4. **Optimize for speed** whenever possible without changing the result.
5. If a PR's branch has already been merged, restart the branch from current `main` before adding new commits; never stack new work on merged history.
6. This session environment has **no cluster access, no real data, and no R install**. Anything written is reasoning-based, not executed. Always say so rather than implying something was tested.
7. Ask on genuine judgment calls (methodology, ambiguous design, anything that changes semantics) instead of guessing. When flagging a bug, state confidence; apparent bugs have previously turned out to compensate for real upstream data quirks (see sister repo history).
8. **Audit downstream effects of every change.** Before opening a PR, trace every object/column the change touches through the rest of the script (names, column order, classes, positional indexing) and fix knock-on breakage in the same PR. Added after PR 2, where the user had to clean up fallout by hand (see Troubleshooting history).

## What this project is

PhenoCycler-Fusion (Akoya) multiplexed immunofluorescence of ALS and control **motor cortex (`mcx`)** and **spinal cord (`sc`)**, 70 channels (69 antibodies + DAPI; see `antibody_targets`). Images are analyzed in **HALO / HaloAI**:

- HaloAI segments/classifies major cell types: neurons (several classes), astrocytes, oligodendrocytes, myeloid cells, T cells; plus vasculature, pTDP-43 inclusions, and C9orf72 DPR (poly-GA / poly-GP) inclusions.
- HALO **Area Quantification FL** gives the mean intensity of every channel within each cell/feature mask → this becomes the single-cell "expression" matrix.
- HALO annotation layers (exported as **geoJSON**) mark gray vs white matter and other anatomical/pathological regions, used to label each object, analogous to how Visium spots were annotated in the sister repo [`jakesboles/als_cns_visium`](https://github.com/jakesboles/als_cns_visium) (`scripts/03b_make_halo_gdfs.py`, `03c_make_halo_feature_gdfs.py`, `04_spot_annotation.R`).

Downstream goals: spatially aware analyses, especially **cellular neighborhoods** and how they shift in ALS (sALS, C9orf72-ALS) vs Control, by tissue and compartment.

Phase 1 (current): import HALO object data, annotate objects with spatial context (GM/WM, proximity to/containment of pathology), and build one single-cell spatial object. Strong preference for **R**.

## Repo contents (as of first commit)

- `sample_key.csv`: slide/section key. `code` = section id (`<run>-<slot>`, e.g. `160-1`), `sample_id` = donor (hyphen-free, e.g. `GBB1607`), `block`, `sn` (section number?), `tissue` (`mcx`/`sc`), `date` (staining run date), unnamed trailing notes column.
- `target_als_demographics_compiled.csv`: donor demographics (same file as in the Visium repo). `Case Number` has hyphens (`GWF18-31`, `GBB-16-07`, `AU-066`); strip hyphens to join on `sample_id`.
- `antibody_targets`: one target per line, 70 lines, no header.
- `scripts/myeloid_prelim.R`, `scripts/pathology_prelim.R`: prototype import → `SpatialExperiment` (myeloid) / `SingleCellExperiment` (pathology) → PCA/UMAP/Harmony → presto DE.
- `scripts/pathology_summary.R`: % area of pTDP-43 / DPR from HALO summary exports, Wilcoxon by group.
- Data dirs `myeloid_quant/`, `pathology_quant/` are gitignored; real data lives on Quest at `/projects/b1169/boles/als_cns_pcf`.

## Cohort observations

- 49 sections, 25 donors. All `sample_id`s match the demographics file after hyphen stripping. Six donors in demographics are not in this cohort: AU-066/072/073/085/087, GBB-23-19.
- Group derivation (shared with Visium repo): `C9orf72 mutation == "Y"` → C9orf72-ALS; else `Clinical Diagnosis == "Control"` → Control; else sALS. Note `ALS-FTD` diagnoses fall into sALS unless C9+.
- Section counts by tissue × group (all sections, incl. failed): mcx Control 10 / sALS 9 / C9 5; sc Control 11 / sALS 8 / C9 6.
- **Failed sections to exclude**: `162-6` (GBB1708 sc, "badly folded") and `162-8` (GBB2311 sc, "broken"). Redone as `170-2` and `168-2`. These are the only donor × tissue duplicates. Exclusion should key off the notes column or an explicit list, not be silently handled by joins.
- 22 of 25 donors have both tissues; GBB1819 (mcx only), GBB2108 (mcx only), GBB1807 (sc only).
- Staining runs (`code` prefix 160–170) are a likely batch variable; `date` is 1:1 with run.

## Data-format observations / gotchas

- **`sample_key.csv` has a UTF-8 BOM.** Depending on platform/locale, `read.csv()` may name the first column `X.U.FEFF.code` or `ï..code` instead of `code`. Safe: `read.csv(..., fileEncoding = "UTF-8-BOM")` or `readr::read_csv()`. If `key[, "code"]` ever errors with "undefined columns", this is why.
- The header ends in a trailing comma → unnamed notes column (`X` in `read.csv`).
- HALO object export filenames start with `JSB<code>_...`; scripts recover `code` by splitting on `_` and dropping `JSB`.
- HALO column naming: per-channel `<Marker>.Average.Positive.Intensity`; area column `Region.Area..μm..` (contains `μ`; encoding-fragile, consider matching by regex). Centroid currently = bounding-box midpoint `(XMin+XMax)/2`, `(YMin+YMax)/2`; units presumably pixels (confirm; PCF 20x is ~0.5 µm/px).
- Marker names after the scripts' cleanup (`HLADR`, `Iba1`, `Mac2Galectin3`, `C1QA`, `H2AX`, `HLAA`, …) don't match `antibody_targets` spelling (`HLA-DR`, `Iba-1`, `Galectin3`, `C1Qa`, `H2A.X`, `HLA-A`). A small explicit mapping table will be needed for a canonical marker list.
- **Selective column loading** (`myeloid_prelim.R`, PR 2): HALO exports are read with `data.table::fread(select = ...)` instead of `read.csv()`. The header is read first (`nrows = 0`) and passed through `make.names(unique = T)` so column names are byte-identical to what `read.csv()` produced; downstream code (`Region.Area..μm..`, `*.Average.Positive.Intensity`, the marker-name cleanup) is unchanged. Kept: `Object.Id`, `XMin/XMax/YMin/YMax`, `Region.Area*`, all `*Average.Positive.Intensity`. `fread(select = )` returns columns in **file order**, so nothing downstream may rename columns positionally. `data.table = F` is required because later code sets `rownames()` on a frame derived from `tab`, which data.tables don't support. Called as `data.table::fread` (no `library(data.table)`) to avoid masking dplyr/purrr (`first`, `last`, `between`, `transpose`). Extend to the other scripts once validated on real data.
- Prototype loops take `feature_cols` from the **last** file read; this assumes every export has identical channel columns in the same order. Worth asserting.

## Preprocessing of intensity data

### Current state (`myeloid_prelim.R`, PR 4)
- Area filter: objects with area ≤ 65 µm² dropped.
- `counts` assay = raw HALO mean intensities (the name is a SCE convention; these are not counts). `asinh` assay = `asinh(counts / cofactor)`, `cofactor <- 5` as a placeholder.
- ComBat-seq and `logNormCounts` removed (with `library(singleCellTK)`/`library(scuttle)`, whose only uses they were). All downstream `combat_*` references now point at `asinh`. `dittoDimPlot(sce, "CD163")` and the TMEM119 `dittoPlot` previously had no `assay =` and so plotted the dittoSeq default assay (likely raw `counts`, despite the "Normalized intensity" label); they now pass `assay = "asinh"` explicitly.
- PCA on `asinh` with `scale = T` (per-marker z-scoring, cf. Hickey 2021), Harmony on `code` for the embedding only.
- presto `logFC` on the `asinh` assay is a difference of mean asinh values (≈ log-ratio for values well above the cofactor).

### Rationale
- **Why asinh**: intensity distributions are heavily right-skewed. asinh(x/c) is ~linear for x ≪ c and ~log(2x/c) for x ≫ c, so it compresses the bright tail like a log while staying defined and stable at 0/near-background values. Standard in cytometry (Bendall 2011; Nowicka 2017 CyTOF workflow, cofactor 5) and multiplexed imaging (Windhager 2023 IMC protocol, cofactor 1 for IMC counts).
- **Cofactor is scale-dependent and must be tuned**: c should sit near the background/noise level of each channel so noise is compressed and positive signal is log-scaled. Cytometry's 5 (CyTOF counts) / ~150 (flow) don't transfer to HALO intensities directly. Plan: inspect per-marker raw distributions (e.g. `apply(counts, 1, quantile, c(.5, .9, .99))` or density plots on log scale) and pick a global or per-marker cofactor. Per-marker works as-is: a vector ordered by `rownames(sce)` divides row-wise under R's column-major recycling.
- **Why drop ComBat-seq**: (1) it models negative-binomial counts (Zhang 2020) and these are continuous intensities; (2) more importantly, it was run with `batch = code`, and every section is one donor in one group, so removing per-section mean differences removes donor/group differences, the signal of interest. Batch correction with batch confounded or unbalanced with group biases downstream group tests (Nygaard 2016).
- **Why logNormCounts is inappropriate**: library-size factors assume each cell's total is technical (sequencing depth). Summed intensity across 70 unrelated antibodies isn't a depth proxy; dividing by it injects cell-type composition into every marker.
- **Batch handling going forward**: correct only the embedding (Harmony on PCA, by `code` or staining run) for clustering/neighborhoods; keep expression values uncorrected for group comparisons and model batch/donor in the statistical test instead. Slide-mean scaling (Harris 2022, mxnorm) was best in their evaluation, but in this design one slide = one donor section, so it has the same confounding problem as ComBat-by-`code`.

### Proposed full pipeline (not yet implemented; for discussion)
1. QC: drop failed sections (`162-6`, `162-8`); area filter; DAPI-low / extreme-area objects (segmentation artifacts); optionally clip per-marker at the 99.9th percentile per section to tame hot pixels/debris.
2. Transform: `asinh(x / cofactor)` with tuned (likely per-marker) cofactors.
3. Embedding/clustering: per-marker z-score (`scale = T`) on a lineage panel → PCA → Harmony (`code`) → clustering/UMAP; cell types otherwise come from HaloAI classes.
4. Spatial: neighborhoods via `imcRtools` with `img_id = "code"`.
5. Group comparisons: never treat cells as replicates (Squair 2021; Zimmerman 2021). Either aggregate per section × cell type (mean asinh) and test with limma/lm (`~ group + sex + age`), or fit mixed models on cells with a donor random effect (`(1 | sample_id)`, plus `code` once both tissues are in). Same logic as diffcyt/CyTOF workflow (Nowicka 2017).

### Open question for user
- HALO's `Average Positive Intensity` may be the mean over **positive (above-threshold) pixels** only, not the mean over the whole object mask. If so, values are conditional on positivity, cells with no positive pixels likely read 0, and thresholds set in HALO shape the data. Confirm which HALO output column is the whole-mask mean, if any.

### References (verified on PubMed)
- Bendall SC et al. 2011 Science. doi:10.1126/science.1198704
- Nowicka M et al. 2017 F1000Research, CyTOF workflow. doi:10.12688/f1000research.11622.3
- Windhager J et al. 2023 Nat Protoc, end-to-end multiplexed imaging workflow (imcRtools). doi:10.1038/s41596-023-00881-0
- Hickey JW et al. 2021 Front Immunol, CODEX normalization vs cell-type accuracy. doi:10.3389/fimmu.2021.727626
- Harris CR et al. 2022 Bioinformatics, slide-to-slide variation in MxIF. doi:10.1093/bioinformatics/btab877
- Zhang Y et al. 2020 NAR Genom Bioinform, ComBat-seq. doi:10.1093/nargab/lqaa078
- Nygaard V et al. 2016 Biostatistics, batch correction with unbalanced groups. doi:10.1093/biostatistics/kxv027
- Squair JW et al. 2021 Nat Commun, pseudoreplication in single-cell DE. doi:10.1038/s41467-021-25960-2
- Zimmerman KD et al. 2021 Nat Commun, mixed models for single-cell pseudoreplication. doi:10.1038/s41467-021-21038-1

## Proposed toolchain (pending user agreement)

- **Object**: `SpatialExperiment` (adopted in `myeloid_prelim.R`, PR 3). Built with `spatialCoordsNames = c("x", "y")`, which moves the centroid columns out of `colData` into `spatialCoords(spe)` (a matrix, columns `x`, `y`). All SCE-based tools in use (scuttle, scater, singleCellTK, harmony, presto, dittoSeq) accept it as an SCE subclass. No images will ever be added (user decision). **`sample_id` vs `code` (decided):** SPE treats `colData$sample_id` as the image/section identifier, but here `sample_id` is the **donor**. This is fine for now because the preliminary data contain one tissue, so each `sample_id` is one section. Once both tissues are loaded, a donor's mcx and sc sections will share a `sample_id` with overlapping coordinates. **User decision: spatially aware tools are always pointed at `code`** (e.g. `imcRtools::buildSpatialGraph(img_id = "code")`); `sample_id` stays the donor id and is not renamed.
- **Import speed**: `data.table::fread` (or `vroom`) with `select=` on only the needed columns; HALO object exports can be very large.
- **geoJSON annotations**: `sf` in R (`st_read` → `st_make_valid` → `st_union` per category; point-in-polygon of cell centroids via `st_intersects`/`st_within`, spatially indexed). Keeps the pipeline in R; avoids the Visium repo's separate Python/geopandas step. HALO exports drawn regions as LineStrings in the Visium project; `st_polygonize`/`st_cast` needed similarly.
- **Pathology proximity**: pTDP-43 / DPR inclusions are objects themselves here, so distance-to-nearest-inclusion via kNN (`RANN`/`dbscan`) is faster than polygon distance.
- **Neighborhoods / spatial stats**: `imcRtools` (`buildSpatialGraph`, `aggregateNeighbors`, `testInteractions`), `lisaClust`, `Banksy`, `spicyR` for group-level spatial comparisons with donor as random effect.

## Troubleshooting history

- **PR 2 fallout (fread change, `myeloid_prelim.R`)**: the original code assigned `colnames(tab) <- c("object.id", "area", "x", "y", "code", feature_col_names)` **positionally**, assuming the old `[, c(other_cols, feature_cols)]` reordering. After switching to `fread(select = )` the extra `XMin/XMax/YMin/YMax` columns and file-order columns broke that mapping. User fixed in `0e0abc3`: dropped the reorder/positional rename, applied the marker-name cleanup to all colnames, and removed `XMin..YMax` after computing centroids. That fix left the ID column named `object_id` (from the in-loop `rename`) while the SCE block still referenced `object.id`; corrected in PR 3. Lesson → rule 8.
- **`library(SpatialExperiment)` fails on the Quest analytic node** (PR 3 follow-up, user commit `e64d98e`). Chain of errors, all from `magick` (imported by SpatialExperiment): (1) `libMagick++-7.Q16HDRI.so.5` not found — `magick` was built against `ImageMagick/7.1.0-31`, a module that can't be loaded on the analytic node; (2) after an rpath rebuild, `libpng15.so.15` not found; (3) then `undefined symbol: xmlNanoHTTPMethod` in `libMagickCore`. Fixes:
  - Rebuilt `magick` on a login node with the module loaded and `LDFLAGS += -Wl,--disable-new-dtags,-rpath,/software/ImageMagick/7.1.0-31/lib,-rpath,/hpc/software/spack_v20d1/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/libpng-1.5.30-txgtu3ltmhsnwbi6xrjmtybycdopf6vh/lib` in `~/.R/Makevars` (line removed afterwards). Diagnose with `ldd <pkg>/libs/<pkg>.so | grep "not found"` run **on the analytic node**.
  - The libxml2 error can't be fixed by rpath: R on the analytic node already has `libxml2` 2.13.4 (2025 Spack stack, no nanohttp) in memory at startup (`grep libxml2 /proc/<pid>/maps`). Fix is the first line of the script: `dyn.load(".../libxml2-2.9.10-.../lib/libxml2.so.2", local = F)`, the in-session equivalent of the user's old `LD_PRELOAD`. Must run before any `library()` call in a fresh session; every entry point that loads SpatialExperiment/magick needs it.
  - No ImageMagick on the 2025 stack was installed properly (only a build tree at `/gpfs/software/2025/ImageMagick/ImageMagick` and someone's conda env pass the `xmlNano` check). Permanent fix would be Quest installing ImageMagick on the 2025 stack; user chose to keep `dyn.load()`.
  - `LD_PRELOAD`/`LD_LIBRARY_PATH` can't be set from inside R (`Sys.setenv`) — the loader reads them at process start.
  - The user's Lmod cache is corrupt (`luac ... spiderT ... unexpected symbol`); `rm -rf ~/.cache/lmod` clears it.
