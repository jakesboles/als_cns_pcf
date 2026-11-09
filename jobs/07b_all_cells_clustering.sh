#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07b_all_cells_clustering
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 240G
#SBATCH --time 24:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%j.log

# dataset-wide z (default): sbatch jobs/07b_all_cells_clustering.sh
# log1p + per-section z: sbatch --export=ALL,PCF_NORM=log_section_z jobs/07b_all_cells_clustering.sh

cd /gpfs/projects/b1169/boles/als_cns_pcf

module purge
module load R/4.5.1

Rscript scripts/07b_all_cells_clustering.R
