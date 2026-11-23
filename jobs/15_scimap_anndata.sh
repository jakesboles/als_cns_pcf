#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 15_scimap_anndata
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 64G
#SBATCH --time 1:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%j.log

# sbatch jobs/15_scimap_anndata.sh
# positions and labels only: --export=ALL,PCF_NO_INTENSITIES=1; also keep raw means: PCF_RAW_LAYER=1

cd /gpfs/projects/b1169/boles/als_cns_pcf

ENV=${PCF_SCIMAP_ENV:-$HOME/.conda/envs/scimap}
PY="$ENV/bin/python"

module purge
"$PY" -c "import sys, scimap, anndata; print(sys.executable, anndata.__version__)" || \
  { echo "No scimap env at $ENV; see envs/scimap.yml"; exit 1; }

"$PY" scripts/15_scimap_anndata.py --workers "$SLURM_CPUS_PER_TASK" \
  ${PCF_NO_INTENSITIES:+--no_intensities} ${PCF_RAW_LAYER:+--raw_layer}
