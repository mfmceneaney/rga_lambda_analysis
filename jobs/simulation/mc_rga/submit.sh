#!/bin/bash

#SBATCH --job-name=RGA_LAMBDA_ANALYSIS_mc_rga
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH -c 1
#SBATCH --mem-per-cpu=2000
##SBATCH --gres=disk:5000
#SBATCH --time=24:00:00

$RGA_LAMBDA_ANALYSIS_HOME/jobs/simulation/mc_rga_sss/job.sh
