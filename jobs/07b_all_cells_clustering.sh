#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07b_all_cells_clustering
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 300G
#SBATCH --time 48:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%j.log

# dataset-wide z (default): sbatch jobs/07b_all_cells_clustering.sh
# log1p + per-section z: sbatch --export=ALL,PCF_NORM=log_section_z jobs/07b_all_cells_clustering.sh
# graph: kNN k = 150 by default (outputs <norm>_knn150/); other k: PCF_K=100; previous SNN k = 15: PCF_GRAPH=snn (outputs <norm>/)

cd /gpfs/projects/b1169/boles/als_cns_pcf

module purge
module load R/4.5.1

Rscript scripts/07b_all_cells_clustering.R
