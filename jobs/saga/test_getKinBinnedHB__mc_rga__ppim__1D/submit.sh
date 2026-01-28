#!/bin/bash

#SBATCH --job-name=saga_getKinBinnedHB
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH -c 1
#SBATCH --mem-per-cpu=4G
##SBATCH --gres=disk:5000
#SBATCH --time=24:00:00

export OUTDIR="$RGA_LAMBDA_ANALYSIS_HOME/jobs/saga/test_getKinBinnedHB__mc_rga__ppim__1D"
export YAML="args.yaml"

echo $OUTDIR
echo $YAML

cd $OUTDIR
ls -lrth
pwd
RGA_LAMBDA_ANALYSIS_SAGA_COMMAND "getKinBinnedHB $YAML"
echo DONE
