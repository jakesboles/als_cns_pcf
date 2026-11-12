#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 11_myeloid_dimred
#SBATCH --array 1-7
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 128G
#SBATCH --time 6:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%A_%a.log

# tasks: 1 asinh, 2 dataset_z, 3 section_z, 4 log_dataset_z, 5 log_section_z, 6 asinh_section_z, 7 ad_combat (order = methods in 11)
# all methods: sbatch jobs/11_myeloid_dimred.sh; some: sbatch --array=1,5 jobs/11_myeloid_dimred.sh
# all channels instead of the panel: sbatch --export=ALL,PCF_FEATURES=all jobs/11_myeloid_dimred.sh (outputs <method>_all_features/)
# number of PCs: PCF_NPCS; supervised UMAP label weight (0-1, default 0.5): PCF_SUP_WEIGHT

cd /gpfs/projects/b1169/boles/als_cns_pcf

module purge
module load R/4.5.1

Rscript scripts/11_myeloid_dimred.R
