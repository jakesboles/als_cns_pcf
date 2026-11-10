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
# graph: SNN k = 15 by default (outputs <norm>/); kNN k = 150: PCF_GRAPH=knn (outputs <norm>_knn150/); other k: PCF_K
# all channels in the PCA instead of the 11-marker panel: PCF_FEATURES=all (outputs <norm>_all_features/; 12 PCs, PCF_NPCS to change)

cd /gpfs/projects/b1169/boles/als_cns_pcf

module purge
module load R/4.5.1

Rscript scripts/07b_all_cells_clustering.R
