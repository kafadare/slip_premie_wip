#!/bin/bash
#SBATCH --job-name=rscript_submit
#SBATCH --output=/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/logs/rscript_submit_%x_%j.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/logs/rscript_submit_%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

module load R

# Rscript file name required

# First argument is the R script
R_SCRIPT=$1
shift  # remove first argument so $@ contains only the additional arguments

echo "Running R script: $R_SCRIPT"
echo "Additional arguments: $@"

Rscript "$R_SCRIPT" $@

#INPUT_SCRIPT=$1

#Rscript "$INPUT_SCRIPT"
