#!/bin/bash
#SBATCH --job-name=render_rmd
#SBATCH --output=/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/logs/render_rmd_%x_%j.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/logs/render_rmd_%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

module load R

# Rmd file name required

INPUT_RMD=$1
OUTPUT_DIR=$2

#if no arg supplied, then uses defaults defined inside render.R and input script
PARAMS_CENTILE=${3:-""}
OUTPUT_BASE=${4:-""}

#conda init
#conda activate slip_premie

# --- Build command depending on which args are supplied ---
if [ -n "$PARAMS_CENTILE" ] && [ -n "$OUTPUT_BASE" ]; then
Rscript render.R "$INPUT_RMD" "$OUTPUT_DIR" "$PARAMS_CENTILE" "$OUTPUT_BASE"
elif [ -n "$PARAMS_CENTILE" ]; then
Rscript render.R "$INPUT_RMD" "$OUTPUT_DIR" "$PARAMS_CENTILE"
else
  Rscript render.R "$INPUT_RMD" "$OUTPUT_DIR"
fi
