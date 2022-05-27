#!/bin/bash
#BSUB -J snakemake
#BSUB -n 1
#BSUB -R rusage[mem=4]
#BSUB -W 24:00
#BSUB -eo logs/cluster/%J.stderr

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isablapi
module load singularity/3.6.2

mkdir logs/
mkdir logs/cluster

snakemake \
  --profile lsf \
  --keep-going \
  --restart-times 1 \
  --singularity-args "--bind /juno/work" \
  --rerun-incomplete \
  --conda-frontend conda