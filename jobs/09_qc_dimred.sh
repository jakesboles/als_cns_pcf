#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 09_qc_dimred
#SBATCH --array 1-6
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 128G
#SBATCH --time 12:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%A_%a.log

# tasks: 1 asinh, 2 dataset_z, 3 section_z, 4 log_dataset_z, 5 log_section_z, 6 asinh_section_z (order = methods in 09)
# all methods: sbatch jobs/09_qc_dimred.sh; some: sbatch --array=1,5 jobs/09_qc_dimred.sh
# all channels instead of the panel: sbatch --export=ALL,PCF_FEATURES=all jobs/09_qc_dimred.sh (outputs <method>_all_features/)
# number of PCs: PCF_NPCS

cd /gpfs/projects/b1169/boles/als_cns_pcf

module purge
module load R/4.5.1

Rscript scripts/09_qc_dimred.R
