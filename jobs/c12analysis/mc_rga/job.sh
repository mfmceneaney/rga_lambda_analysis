#!/bin/bash                                                  

export INFILE="/path/test.hipo"
export OUTDIR="$RGA_LAMBDA_ANALYSIS_DIR/jobs/c12analysis/mc_rga"
export name=`echo $INFILE | xargs -n 1 basename`
mkdir -p $OUTDIR
cd $OUTDIR

$C12ANALYSIS/bin/run.sh $INFILE -ch 2212,-211 -be $BEAM_ENERGY_RGA -tpid $TPID_RGA -rn -en -ang -vtx -lk -ik -ma -f -out $OUTDIR/skim_${name}.root -momc Fall2018 1 1

echo DONE