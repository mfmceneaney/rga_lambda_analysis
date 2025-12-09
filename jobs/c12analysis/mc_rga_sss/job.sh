#!/bin/bash                                                  

export INFILE='/path/test.hipo'
export OUTDIR="$RGA_LAMBDA_ANALYSIS_VOL_DIR/jobs/c12analysis/mc_rga_sss"
export name=`echo $INFILE | xargs -n 1 basename`
mkdir -p $OUTDIR
cd $OUTDIR

# Get proton pion skim for lambda
$RGA_LAMBDA_ANALYSIS_C12ANALYSIS_COMMAND $INFILE -ch 2212,-211 -be $RGA_LAMBDA_ANALYSIS_BEAM_ENERGY_RGA -tpid $RGA_LAMBDA_ANALYSIS_TPID_RGA -rn -en -ang -vtx -ik -lk -ma -f -out $OUTDIR/skim_ppim_${name}.root

echo DONE
