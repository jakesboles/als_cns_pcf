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

# submit: sbatch jobs/07_embedding.sh  (other R module: sbatch --export=ALL,R_MODULE=R/4.5.x jobs/07_embedding.sh)

cd /gpfs/projects/b1169/boles/als_cns_pcf

module purge
module load ${R_MODULE:-R/4.5.0} || { echo "Could not load ${R_MODULE:-R/4.5.0}; check 'module avail R'"; exit 1; }
Rscript -e 'cat(R.version.string, .libPaths(), sep = "\n")'

Rscript scripts/07_embedding.R
