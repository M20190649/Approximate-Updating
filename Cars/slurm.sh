#!/bin/env/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=NTCarFC
#SBATCH --array=1-873

module load R/3.4.1
module load boost/1.6.3

R CMD BATCH slurmFC.R ./output/out${SLURM_ARRAY_JOB_ID}.Rout




