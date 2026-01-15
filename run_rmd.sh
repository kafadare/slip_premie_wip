#!/bin/bash
#SBATCH --job-name=render_rmd
#SBATCH --output=render_rmd_%j.out
#SBATCH --error=render_rmd_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

module load R

# Rmd file name required

INPUT_RMD=$1

#if no arg supplied, then uses defaults defined inside render.R
OUTPUT_DIR=${2:-""}
OUTPUT_BASE=${3:-""}

#conda init
#conda activate slip_premie

# --- Build command depending on which args are supplied ---
if [ -n "$OUTPUT_DIR" ] && [ -n "$OUTPUT_BASE" ]; then
Rscript render.R "$INPUT_RMD" "$OUTPUT_DIR" "$OUTPUT_BASE"
elif [ -n "$OUTPUT_DIR" ]; then
Rscript render.R "$INPUT_RMD" "$OUTPUT_DIR"
else
  Rscript render.R "$INPUT_RMD"
fi
