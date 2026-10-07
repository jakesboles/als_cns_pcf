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
- `scripts/myeloid_prelim.R`, `scripts/pathology_prelim.R`: prototype import → `SingleCellExperiment` → PCA/UMAP/Harmony → presto DE.
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
- **Selective column loading** (`myeloid_prelim.R`, PR 2): HALO exports are read with `data.table::fread(select = ...)` instead of `read.csv()`. The header is read first (`nrows = 0`) and passed through `make.names(unique = T)` so column names are byte-identical to what `read.csv()` produced; downstream code (`Region.Area..μm..`, `*.Average.Positive.Intensity`, the marker-name cleanup) is unchanged. Kept: `Object.Id`, `XMin/XMax/YMin/YMax`, `Region.Area*`, all `*Average.Positive.Intensity`. `data.table = F` is required because later code sets `rownames()` on a frame derived from `tab`, which data.tables don't support. Called as `data.table::fread` (no `library(data.table)`) to avoid masking dplyr/purrr (`first`, `last`, `between`, `transpose`). Extend to the other scripts once validated on real data.
- Prototype loops take `feature_cols` from the **last** file read; this assumes every export has identical channel columns in the same order. Worth asserting.

## Methodology notes from prototype scripts (for discussion, not yet acted on)

- Myeloid prototype drops objects with area ≤ 65 µm².
- `runComBatSeq` (count-based, negative binomial) and `logNormCounts` (library-size normalization) are applied to mean intensities. Both assume count data; for MxIF intensities, common choices are arcsinh/log1p transform with per-marker scaling, and Harmony (or nothing) for batch. Raise with the user before changing.
- Harmony on `code` (section) is the current batch strategy. In the Visium repo the user found CCA beat Harmony; unknown whether that will hold here.

## Proposed toolchain (pending user agreement)

- **Object**: `SpatialExperiment` (extends the `SingleCellExperiment` already used; centroids in `spatialCoords`, one object for all sections, `sample_id` = section `code`).
- **Import speed**: `data.table::fread` (or `vroom`) with `select=` on only the needed columns; HALO object exports can be very large.
- **geoJSON annotations**: `sf` in R (`st_read` → `st_make_valid` → `st_union` per category; point-in-polygon of cell centroids via `st_intersects`/`st_within`, spatially indexed). Keeps the pipeline in R; avoids the Visium repo's separate Python/geopandas step. HALO exports drawn regions as LineStrings in the Visium project; `st_polygonize`/`st_cast` needed similarly.
- **Pathology proximity**: pTDP-43 / DPR inclusions are objects themselves here, so distance-to-nearest-inclusion via kNN (`RANN`/`dbscan`) is faster than polygon distance.
- **Neighborhoods / spatial stats**: `imcRtools` (`buildSpatialGraph`, `aggregateNeighbors`, `testInteractions`), `lisaClust`, `Banksy`, `spicyR` for group-level spatial comparisons with donor as random effect.

## Troubleshooting history

(none yet)
