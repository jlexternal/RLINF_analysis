#!/bin/bash
#SBATCH --job-name=RLINF_pilot7
#SBATCH --output=./out/RLINF_pilot7%A_%a.out
#SBATCH --error=./err/RLINF_pilot7%A_%a.err
#SBATCH --partition=fastgen,dellgen,firstgen,gnt,lastgen,secondgen
#SBATCH --array=1-23

cd /home/jlee/RLINF
module load MATLAB/R2018b
matlab -nodesktop -r slurm_run_fit_noisyINF_rlinf
