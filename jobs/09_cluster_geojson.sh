#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 09_cluster_geojson
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 64G
#SBATCH --time 2:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_pcf/logs/%x_%j.log

# per-type clusters from 07 (astrocyte, cortical_neuron, myeloid, oligodendrocyte, tcell):
#   sbatch --export=ALL,PCF_GROUP=oligodendrocyte jobs/09_cluster_geojson.sh
# without PCF_GROUP: 08's all-cell clusters
# all images, all clusters at Leiden 0.3: sbatch jobs/09_cluster_geojson.sh
# other resolution / some clusters / some images:
#   sbatch --export=ALL,PCF_RES=0.5,PCF_CLUSTERS="3 7",PCF_CODE="160-1 160-2" jobs/09_cluster_geojson.sh

cd /gpfs/projects/b1169/boles/als_cns_pcf

ENV=${PCF_ENV:-$HOME/.conda/envs/pcf_quant}

PY="$ENV/bin/python"

module purge
module load mamba/23.1.0
eval "$(conda shell.bash hook)"
conda activate "$ENV" || { echo "Could not activate $ENV"; exit 1; }
"$PY" -c "import sys, numpy; print(sys.executable, numpy.__version__)" || exit 1

"$PY" scripts/09_cluster_geojson.py --workers "$SLURM_CPUS_PER_TASK" --resolution "${PCF_RES:-0.3}" ${PCF_GROUP:+--group $PCF_GROUP} ${PCF_CLUSTERS:+--clusters $PCF_CLUSTERS} ${PCF_CODE:+--codes $PCF_CODE}
