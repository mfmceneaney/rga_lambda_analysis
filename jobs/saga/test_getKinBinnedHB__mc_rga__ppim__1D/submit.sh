#!/bin/bash

#SBATCH --job-name=mc_rga_saga_getKinBinnedAsym
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=gpu
#SBATCH --account=clas12
#SBATCH -c 8
#SBATCH --mem-per-cpu=2G
#SBATCH --gres=disk:5000
#SBATCH --time=24:00:00
#SBATCH --mail-user=matthew.mceneaney@duke.edu

export MYEXECUTABLE=$SAGA_BUILD_DIR/saga/getKinBinnedAsym
export OUTDIR=$RGA_LAMBDA_ANALYSIS/jobs/saga/test_getKinBinnedAsym__mc_rga__ppim__1D/
export YAML=args.yaml

echo $MYEXECUTABLE
echo $OUTDIR
echo $YAML

cd $OUTDIR
ls -lrth
pwd
$MYEXECUTABLE $YAML
echo DONE
