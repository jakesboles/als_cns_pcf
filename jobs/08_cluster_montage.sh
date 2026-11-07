#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 08_cluster_montage
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 64G
#SBATCH --time 2:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%j.log

# all 07 groups at Leiden 0.3: sbatch jobs/08_cluster_montage.sh
# options: sbatch --export=ALL,PCF_RES=0.5,PCF_GROUP="myeloid oligodendrocyte",PCF_N=50 jobs/08_cluster_montage.sh

cd /gpfs/projects/b1169/boles/als_cns_pcf

ENV=${PCF_ENV:-$HOME/.conda/envs/pcf_quant}

PY="$ENV/bin/python"

module purge
module load mamba/23.1.0
eval "$(conda shell.bash hook)"
conda activate "$ENV" || { echo "Could not activate $ENV"; exit 1; }
"$PY" -c "import sys, numpy, PIL; print(sys.executable, numpy.__version__, PIL.__version__)" || exit 1

"$PY" scripts/08_cluster_montage.py --workers "$SLURM_CPUS_PER_TASK" --resolution "${PCF_RES:-0.3}" ${PCF_GROUP:+--groups $PCF_GROUP} ${PCF_N:+--n_cells $PCF_N}
