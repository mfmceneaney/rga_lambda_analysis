#!/bin/bash

## IMPORTANT! Load local CCDB copy.  Make sure these are in your startup script.
# module load sqlite/dev
# module load xrootd

# Set job info
export MCINDEX=0 #NOTE: This will be set by `setup.sh`.
export INFILE="file.root"
export TREE=$RGA_LAMBDA_ANALYSIS_TREE
export ROOT_SCRIPT="$RGA_LAMBDA_ANALYSIS_HOME/macros/TreeToLund.C"
export PREFIX="out_"
export BASENAME="${PREFIX}${MCINDEX}"
export NEVENTS=$RGA_LAMBDA_ANALYSIS_NMAX
export GCARD=$RGA_LAMBDA_ANALYSIS_GCARD_RGA_GEMC
export YAML=$RGA_LAMBDA_ANALYSIS_YAML_RGA_COAT
export OUTDIR="$RGA_LAMBDA_ANALYSIS_VOL_DIR/jobs/simulation/mc_rga_sss"

# Set parameters for TreeToLund macro
export nTargetNucleons=1
export nTargetProtons=1
export targetPol=0
export beamPol=1
export beamPid=11
export beamEnergy=$RGA_LAMBDA_ANALYSIS_BEAM_ENERGY_RGA
export targetPid=$RGA_LAMBDA_ANALYSIS_TPID_RGA
export processId=1

mkdir -p $OUTDIR

function check_task_status() {

    # Set arguments
    export TASK_NAME=$1
    export TASK_OUTPUT=$2
    export TASK_STATUS=$3
    export TASK_EXIT_CODE=$4

    # Print info
    echo "Task name:   $TASK_NAME"
    echo "Task output: $TASK_OUTPUT"
    echo "Task status: $TASK_STATUS"

    # Exit if status is nonzero
    if [ $TASK_STATUS -ne 0 ]; then
        echo "Task failed. Exit code $TASK_EXIT_CODE."
        exit $TASK_EXIT_CODE
    fi
    
}

# Convert ROOT String Spinner file to lund file
export OUTDIR_LUND="${OUTDIR}/lund"
mkdir -p $OUTDIR_LUND
cd $OUTDIR_LUND
export LUNDFILE=$OUTDIR_LUND/${BASENAME}.lund
export ROOT_COMMAND="${ROOT_SCRIPT}(\\\"${INFILE}\\\",\\\"${TREE}\\\",\\\"${LUNDFILE}\\\", \
$MCINDEX, $NEVENTS, \
$nTargetNucleons, $nTargetProtons, $targetPol, $beamPol, $beamPid, $beamEnergy, $targetPid, $processId)"
echo EXECUTING: "RGA_LAMBDA_ANALYSIS_SAGA_COMMAND \"root -l -b -q \\\"$ROOT_COMMAND\\\"\""
RGA_LAMBDA_ANALYSIS_SAGA_COMMAND "root -l -b -q \"$ROOT_COMMAND\""
export LUND_TASK_STATUS=$?
cd - #NOTE: cd back to wherever you were before
check_task_status "$ROOT_SCRIPT" $LUNDFILE $LUND_TASK_STATUS 1

# Run GEMC detector simulation
export OUTDIR_GEMC="${OUTDIR}/gemc"
mkdir -p $OUTDIR_GEMC
export GEMCFILE=$OUTDIR_GEMC/$BASENAME.hipo
echo EXECUTING: RGA_LAMBDA_ANALYSIS_GEMC_COMMAND \"$GCARD -SAVE_ALL_MOTHERS=1 -SKIPREJECTEDHITS=1 -NGENP=50 -INTEGRATEDRAW=\\\"\*\\\" -USE_GUI=0 -RUNNO=11 -INPUT_GEN_FILE=\\\"LUND, $LUNDFILE\\\" -OUTPUT=\\\"hipo, $GEMCFILE\\\" -N=$NEVENTS\"
RGA_LAMBDA_ANALYSIS_GEMC_COMMAND "$GCARD -SAVE_ALL_MOTHERS=1 -SKIPREJECTEDHITS=1 -NGENP=50 -INTEGRATEDRAW=\"*\" -USE_GUI=0 -RUNNO=11 -INPUT_GEN_FILE=\"LUND, $LUNDFILE\" -OUTPUT=\"hipo, $GEMCFILE\" -N=$NEVENTS"
export GEMC_TASK_STATUS=$?
check_task_status "gemc" $GEMCFILE $GEMC_TASK_STATUS 2

# Run CLAS12 reconstruction
export OUTDIR_REC="${OUTDIR}/cooked"
mkdir -p $OUTDIR_REC
export RECFILE=$OUTDIR_REC/$BASENAME.hipo
echo EXECUTING: RGA_LAMBDA_ANALYSIS_RECON_UTIL_COMMAND \"-i $GEMCFILE -o $RECFILE -n $NEVENTS -y $YAML\"
RGA_LAMBDA_ANALYSIS_RECON_UTIL_COMMAND "-i $GEMCFILE -o $RECFILE -n $NEVENTS -y $YAML"
export REC_TASK_STATUS=$?
check_task_status "recon-util" $RECFILE $REC_TASK_STATUS 3

# Create DST
export OUTDIR_DST="${OUTDIR}/dst"
mkdir -p $OUTDIR_DST
export DSTFILE=$OUTDIR_DST/$BASENAME.hipo
export FILTER_BANKS='RUN::*,RAW::epics,RAW::scaler,HEL::flip,HEL::online,REC::*,RECFT::*,MC::RecMatch,MC::GenMatch,MC::Particle,MC::User,MC::Header,MC::Lund,MC::Event'
$RGA_LAMBDA_ANALYSIS_HIPO_UTILS_COMMAND "-filter -b \"$FILTER_BANKS\" -merge -o $DSTFILE $RECFILE"
export DST_TASK_STATUS=$?
check_task_status "hipo-utils -filter" $DSTFILE $DST_TASK_STATUS 4

# Clean up intermediate files
#rm $LUNDFILE #NOTE: Keep these files because you need to know how much XS you generated and that is recorded in the file names.
rm $GEMCFILE
rm $RECFILE
