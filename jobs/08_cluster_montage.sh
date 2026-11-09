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

# all 07 groups at every Leiden resolution in 07: sbatch jobs/08_cluster_montage.sh
# options: sbatch --export=ALL,PCF_RES="0.3 1",PCF_GROUP="myeloid oligodendrocyte",PCF_N=50 jobs/08_cluster_montage.sh
# 07b (all cells): sbatch --export=ALL,PCF_CLUSTERS_DIR=data/07b_all_cells,PCF_OUT=results/08_cluster_montage_07b/ jobs/08_cluster_montage.sh
# separate .jpg files instead of images embedded in the HTML: add PCF_NO_EMBED=1

cd /gpfs/projects/b1169/boles/als_cns_pcf

ENV=${PCF_ENV:-$HOME/.conda/envs/pcf_quant}

PY="$ENV/bin/python"

module purge
module load mamba/23.1.0
eval "$(conda shell.bash hook)"
conda activate "$ENV" || { echo "Could not activate $ENV"; exit 1; }
"$PY" -c "import sys, numpy, PIL; print(sys.executable, numpy.__version__, PIL.__version__)" || exit 1

"$PY" scripts/08_cluster_montage.py --workers "$SLURM_CPUS_PER_TASK" ${PCF_RES:+--resolutions $PCF_RES} ${PCF_GROUP:+--groups $PCF_GROUP} ${PCF_N:+--n_cells $PCF_N} ${PCF_NO_EMBED:+--no_embed} ${PCF_CLUSTERS_DIR:+--clusters_dir $PCF_CLUSTERS_DIR} ${PCF_OUT:+--results_dir $PCF_OUT}
