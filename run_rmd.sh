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
SAVE_MARKDOWN_FLAG=$3
SAVE_MARKDOWN_NAME=$4

shift 4

echo "Running R script: $R_SCRIPT"
echo "Output Directory: $OUTPUT_DIR"
echo "Save Markdown Flag: $SAVE_MARKDOWN_FLAG"
echo "Save Markdown Name: $SAVE_MARKDOWN_NAME"
echo "Additional arguments: $@"

Rscript render.R "$INPUT_RMD" "$OUTPUT_DIR" $SAVE_MARKDOWN_FLAG $SAVE_MARKDOWN_NAME $@

