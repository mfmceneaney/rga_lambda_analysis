#!/bin/bash

#SBATCH --job-name=saga_getBinKinematics_pi
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH -c 4
#SBATCH --mem-per-cpu=2G
##SBATCH --gres=disk:1000
#SBATCH --time=24:00:00

export OUTDIR="$RGA_LAMBDA_ANALYSIS_HOME/jobs/saga/test_getBinKinematics__mc_rga__ppim"
export YAML="args.yaml"

echo $OUTDIR
echo $YAML

cd $OUTDIR
ls -lrth
pwd
RGA_LAMBDA_ANALYSIS_SAGA_COMMAND "getBinKinematics $YAML"
echo DONE
