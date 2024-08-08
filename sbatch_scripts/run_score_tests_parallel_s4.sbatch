#!/bin/bash
#SBATCH --job-name constraint_inf_parallel_s4
#SBATCH --partition short
#SBATCH --ntasks 1                   
#SBATCH --time 0-11:50                
#SBATCH --mem-per-cpu=5000     
#SBATCH --array=1-758
#SBATCH -o constraint_inf_full_s4_%j.out           
#SBATCH -e constraint_inf_full_s4_%j.err            
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=teichs@uw.edu 
module load R
Rscript constraint_effects/scripts/run_score_tests_parallel_s4.R $SLURM_ARRAY_TASK_ID