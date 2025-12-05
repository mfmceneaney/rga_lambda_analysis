#!/bin/bash

# Modify submit scripts with farm out and partition info
for file in $RGA_LAMBDA_ANALYSIS_HOME/jobs/*/*/*submit.sh; do
    sed -i.bak "s;/farm_out/%u;$RGA_LAMBDA_ANALYSIS_FARM_OUT;g" $file
    sed -i.bak "s;partition=production;partition=$RGA_LAMBDA_ANALYSIS_HPC_PARTITION;g" $file
    if [ -z "$RGA_LAMBDA_ANALYSIS_HPC_ACCOUNT" ]; then
        sed -i.bak "s;#SBATCH --account=;##SBATCH --account=;g" $file
    else
        sed -i.bak "s;#SBATCH --account=.*;#SBATCH --account=$RGA_LAMBDA_ANALYSIS_HPC_ACCOUNT;" $file
    fi
done

# Modify argument yamls
for file in $RGA_LAMBDA_ANALYSIS_HOME/jobs/saga/*/args*.yaml; do
    sed -i.bak "s;RGA_LAMBDA_ANALYSIS_PDG_ALPHA;$RGA_LAMBDA_ANALYSIS_PDG_ALPHA;g" $file
    sed -i.bak "s;DT_RGA_BEAM_POLARIZATION;$DT_RGA_BEAM_POLARIZATION;g" $file
    sed -i.bak "s;MC_RGA_BEAM_POLARIZATION;$MC_RGA_BEAM_POLARIZATION;g" $file
done

# Create output directories
mkdir -p $RGA_LAMBDA_ANALYSIS_FARM_OUT
mkdir -p $RGA_LAMBDA_ANALYSIS_VOL_DIR
