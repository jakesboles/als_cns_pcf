#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 01_quantify_cells
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 120G
#SBATCH --time 4:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%A_%a.log

# submit: sbatch --array=1-$(ls pcf_images/*.qptiff | wc -l)%20 jobs/01_quantify_cells.sh

cd /gpfs/projects/b1169/boles/als_cns_pcf

module purge
module load python-miniconda3/4.12.0
source activate pcf_quant

IMAGE=$(ls pcf_images/*.qptiff | sed -n "${SLURM_ARRAY_TASK_ID}p")

python3 scripts/01_quantify_cells.py --image "$IMAGE" --workers "$SLURM_CPUS_PER_TASK"
