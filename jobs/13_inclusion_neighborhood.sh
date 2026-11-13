#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 13_inclusion_neighborhood
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 2:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%j.log

# all sections in 09 filtered_data.csv: sbatch jobs/13_inclusion_neighborhood.sh
# selected sections: sbatch --export=ALL,PCF_CODE="160-1 160-2" jobs/13_inclusion_neighborhood.sh
# other radii (um): add PCF_RADII="5 15 30 60"

cd /gpfs/projects/b1169/boles/als_cns_pcf

ENV=${PCF_ENV:-$HOME/.conda/envs/pcf_quant}

PY="$ENV/bin/python"

module purge
module load mamba/23.1.0
eval "$(conda shell.bash hook)"
conda activate "$ENV" || { echo "Could not activate $ENV"; exit 1; }
"$PY" -c "import sys, numpy; print(sys.executable, numpy.__version__)" || exit 1

"$PY" scripts/13_inclusion_neighborhood.py --workers "$SLURM_CPUS_PER_TASK" ${PCF_CODE:+--codes $PCF_CODE} ${PCF_RADII:+--radii $PCF_RADII}
