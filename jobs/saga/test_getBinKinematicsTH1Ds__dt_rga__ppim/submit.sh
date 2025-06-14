#!/bin/bash

#SBATCH --job-name=dt_rga_saga_getBinKinematicsTH1Ds
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH -c 4
#SBATCH --mem-per-cpu=2G
#SBATCH --gres=disk:1000
#SBATCH --time=8:00:00

export MYEXECUTABLE=$SAGA_BUILD_DIR/saga/getBinKinematicsTH1Ds
export OUTDIR=$RGA_LAMBDA_ANALYSIS_HOME/jobs/saga/test_getBinKinematicsTH1Ds__dt_rga__ppim
export YAML=args.yaml

mkdir -p $OUTDIR
cd $OUTDIR
$MYEXECUTABLE $YAML
echo DONE
