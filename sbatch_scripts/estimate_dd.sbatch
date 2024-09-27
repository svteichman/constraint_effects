#!/bin/bash
#SBATCH --job-name estimate_dd
#SBATCH --partition medium
#SBATCH --ntasks 1                   
#SBATCH --time 6-11:50                
#SBATCH --mem-per-cpu=5000             
#SBATCH --array=1-84
#SBATCH -o estimate_dd_%j.out           
#SBATCH -e estimate_dd_%j.err            
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=teichs@uw.edu 
module load R
Rscript constraint_effects/scripts/estimate_dd.R $SLURM_ARRAY_TASK_ID