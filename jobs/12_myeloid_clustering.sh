#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 12_myeloid_clustering
#SBATCH --array 1-2
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 128G
#SBATCH --time 6:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%A_%a.log

# tasks: 1 log_section_z, 2 dataset_z (order = runs in 12); needs 11's output for each
# any other 11 run: sbatch --export=ALL,PCF_METHOD=asinh_section_z jobs/12_myeloid_clustering.sh (with --array=1)
# number of Harmony dimensions: PCF_NPCS

cd /gpfs/projects/b1169/boles/als_cns_pcf

module purge
module load R/4.5.1

Rscript scripts/12_myeloid_clustering.R
