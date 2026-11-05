#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07_embedding
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 128G
#SBATCH --time 12:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%j.log

cd /gpfs/projects/b1169/boles/als_cns_pcf

module purge
module load R/4.5.1

Rscript scripts/07_embedding.R
