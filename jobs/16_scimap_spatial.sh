#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 16_scimap_spatial
#SBATCH --array 1-2
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 96G
#SBATCH --time 12:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%A_%a.log

# both regions (task 1 GM, task 2 WM): sbatch jobs/16_scimap_spatial.sh
# one region: sbatch --array=1 jobs/16_scimap_spatial.sh, or --export=ALL,PCF_REGION=WM
# options: PCF_RADII="15 30 60", PCF_NPERM=1000, PCF_K=10, PCF_CODE="160-1 160-2", PCF_NO_NBHD=1

cd /gpfs/projects/b1169/boles/als_cns_pcf

REGIONS=(GM WM)
REGION=${PCF_REGION:-${REGIONS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}}
[ -n "$REGION" ] || { echo "No region for array task $SLURM_ARRAY_TASK_ID"; exit 1; }

ENV=${PCF_SCIMAP_ENV:-$HOME/.conda/envs/scimap}
PY="$ENV/bin/python"

module purge
"$PY" -c "import sys, scimap, anndata; print(sys.executable, anndata.__version__)" || \
  { echo "No scimap env at $ENV; see envs/scimap.yml"; exit 1; }

MPLBACKEND=Agg "$PY" scripts/16_scimap_spatial.py --region "$REGION" \
  ${PCF_RADII:+--radii $PCF_RADII} ${PCF_NPERM:+--n_perm $PCF_NPERM} ${PCF_K:+--k $PCF_K} \
  ${PCF_CODE:+--codes $PCF_CODE} ${PCF_NO_NBHD:+--no_neighborhoods}
