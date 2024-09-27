#!/bin/bash
#SBATCH --job-name score_tests_dd
#SBATCH --partition medium
#SBATCH --ntasks 1                   
#SBATCH --time 6-11:50                
#SBATCH --mem-per-cpu=10000             
#SBATCH --array=1-84
#SBATCH -o score_tests_dd_%j.out           
#SBATCH -e score_tests_dd_%j.err            
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=teichs@uw.edu 
module load R
Rscript constraint_effects/scripts/run_score_tests_dd.R $SLURM_ARRAY_TASK_ID