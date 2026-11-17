#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 14_cell_pair_spatial
#SBATCH --array 1-2
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 32G
#SBATCH --time 6:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%A_%a.log

# both regions (task 1 GM, task 2 WM): sbatch jobs/14_cell_pair_spatial.sh
# one region: sbatch --array=2 jobs/14_cell_pair_spatial.sh, or --export=ALL,PCF_REGION=WM
# options: PCF_CODE="160-1 160-2", PCF_RADII="15 30 60", PCF_NPERM=1000

cd /gpfs/projects/b1169/boles/als_cns_pcf

REGIONS=(GM WM)
REGION=${PCF_REGION:-${REGIONS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}}
[ -n "$REGION" ] || { echo "No region for array task $SLURM_ARRAY_TASK_ID"; exit 1; }

ENV=${PCF_ENV:-$HOME/.conda/envs/pcf_quant}

PY="$ENV/bin/python"

module purge
module load mamba/23.1.0
eval "$(conda shell.bash hook)"
conda activate "$ENV" || { echo "Could not activate $ENV"; exit 1; }
"$PY" -c "import sys, numpy, scipy; print(sys.executable, numpy.__version__, scipy.__version__)" || \
  { echo "scipy missing: mamba install -p $ENV -c conda-forge scipy"; exit 1; }

"$PY" scripts/14_cell_pair_spatial.py --region "$REGION" --workers "$SLURM_CPUS_PER_TASK" \
  ${PCF_CODE:+--codes $PCF_CODE} ${PCF_RADII:+--radii $PCF_RADII} ${PCF_NPERM:+--n_perm $PCF_NPERM}
