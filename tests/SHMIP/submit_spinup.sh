#!/bin/bash -l
#SBATCH --job-name=baseline
#SBATCH --account=def-c3dow    # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=0-01:00         # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # adjust this if you are using parallel commands
#SBATCH --mem=8G               # adjust this according to the memory requirement per node you need
#SBATCH --output=spinup.log

module load StdEnv/2020

# Choose a version of MATLAB by loading a module:
module load matlab/2020b

# Remove -singleCompThread below if you are using parallel commands:
matlab -nodisplay -singleCompThread -r "spinup_baseline; exit"
