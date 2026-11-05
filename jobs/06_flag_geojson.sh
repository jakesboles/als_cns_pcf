#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 06_flag_geojson
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 64G
#SBATCH --time 2:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%j.log

# all images in 05 output: sbatch jobs/06_flag_geojson.sh
# selected images: sbatch --export=ALL,PCF_CODE="160-1 160-2" jobs/06_flag_geojson.sh
# only dropped cells: add PCF_DROPPED_ONLY=1 to --export

cd /gpfs/projects/b1169/boles/als_cns_pcf

ENV=${PCF_ENV:-$HOME/.conda/envs/pcf_quant}

PY="$ENV/bin/python"

module purge
module load mamba/23.1.0
eval "$(conda shell.bash hook)"
conda activate "$ENV" || { echo "Could not activate $ENV"; exit 1; }
"$PY" -c "import sys, numpy; print(sys.executable, numpy.__version__)" || exit 1

"$PY" scripts/06_flag_geojson.py --workers "$SLURM_CPUS_PER_TASK" ${PCF_CODE:+--codes $PCF_CODE} ${PCF_DROPPED_ONLY:+--dropped_only}
