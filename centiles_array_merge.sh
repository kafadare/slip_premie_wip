#!/bin/bash
#SBATCH --job-name=array-job
#SBATCH --output=/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/logs/%x_%j.out
#SBATCH -t 0-6:00

###sample usage
##NO_PHENO -- number of phenotypes/number of array jobs (260 for all aparc, global, and aseg phenotypes)
##ARRAY SCRIPT -- script that runs on each iteration of array (in this case to calculate centiles for each pheno) get_centiles.r, pass full path
##FINAL SCRIPT -- script that runs at the end, in this case merges phenos combine_centile_csv.r, pass full path
##MODEL DIR -- folder for gamlss fits: "/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/RDS/" 
##DATA_FILENAME -- file with raw values: "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/CSV/slip_median.csv"
##OUT_FOLDER -- where to save output centiles: "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/CSV/centile_csv/slip-median-all/"
##MODEL_TYPE -- which type of model (should ideally match the file above) "median"
##JOB_NAME -- any job name desired, defaults to "array-job"

NO_PHENO=$1
ARRAY_SCRIPT=$2
FINAL_SCRIPT=$3
MODEL_DIR=$4
DATA_FILENAME=$5
OUT_FOLDER=$6
MODEL_TYPE=$7
JOB_NAME={$8:-"array-job"}

echo "Arguments: $@"

job_array_id=$(sbatch --array=1-$NO_PHENO -J $JOB_NAME rscript-submit-insidearray.sh $ARRAY_SCRIPT $MODEL_DIR $DATA_FILENAME $OUT_FOLDER $MODEL_TYPE | awk '{print $4}')
echo "Job Array ID: $job_array_id"

sbatch --dependency=afterok:$job_array_id rscript-submit.sh $FINAL_SCRIPT $OUT_FOLDER $MODEL_TYPE