#!/bin/env bash
#SBATCH --job-name=carsFC
#SBATCH --ntasks=1
#SBATCH --array=1-10
#SBATCH --output=output.txt
module load R/3.4.1
module load boost/1.63.0
module load gcc/4.9.1
R --vanilla < slurmFCAll.R

