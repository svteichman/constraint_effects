#!/bin/bash
#SBATCH --job-name constraint_est_mean
#SBATCH --partition short
#SBATCH --ntasks 1                   
#SBATCH --time 0-11:50                
#SBATCH --mem-per-cpu=5000             
#SBATCH -o constraint_est_mean_%j.out           
#SBATCH -e constraint_est_mean_%j.err            
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=teichs@uw.edu 
module load R
Rscript constraint_effects/scripts/estimate_mean.R 