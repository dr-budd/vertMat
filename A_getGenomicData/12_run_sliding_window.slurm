#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=2-00:00:00
#SBATCH --output=%x.%jx.log

## load R
module load R/4.1.3

## move to wd
cd $SLURM_SUBMIT_DIR

## print the start time
echo "R script started running at $(date)"

## run R script
Rscript 12_sliding_window.R

## print the end time
echo "R script finished running at $(date)"

## end script
