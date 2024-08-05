#!/bin/bash
#SBATCH --job-name constraint_inf
#SBATCH --partition medium
#SBATCH --ntasks 1                   
#SBATCH --time 6-11:50                
#SBATCH --mem-per-cpu=10000             
#SBATCH --array=1-45
#SBATCH -o constraint_inf_%j.out           
#SBATCH -e constraint_inf_%j.err            
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=teichs@uw.edu 
module load R
Rscript constraint_effects/scripts/run_score_tests.R $SLURM_ARRAY_TASK_ID