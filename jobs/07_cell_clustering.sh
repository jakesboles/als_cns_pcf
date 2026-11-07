#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07_cell_clustering
#SBATCH --array 1-5
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 128G
#SBATCH --time 12:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%A_%a.log

# tasks: 1 astrocyte, 2 cortical_neuron, 3 myeloid, 4 oligodendrocyte, 5 tcell (cd4 + cd8); order = cell_groups in 07
# all groups: sbatch jobs/07_cell_clustering.sh
# some groups: sbatch --array=3,5 jobs/07_cell_clustering.sh

cd /gpfs/projects/b1169/boles/als_cns_pcf

module purge
module load R/4.5.1

Rscript scripts/07_cell_clustering.R
