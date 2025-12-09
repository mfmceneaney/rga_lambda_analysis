#!/bin/bash

#----- DEFAULT VARIABLES -----#

# Set variables for string spinner simulation
export RGA_LAMBDA_ANALYSIS_TREE="tree" #NOTE: CHANGE AS NEEDED.
export RGA_LAMBDA_ANALYSIS_NEVENTS=10000000 # Number of events in each file #NOTE: CHANGE AS NEEDED.
export RGA_LAMBDA_ANALYSIS_NMAX=10000 # Maximum number of events allowed in each output lund file

# Set beam energies and target lund pids for CLAS12-Analysis #NOTE: CHANGE AS NEEDED
export RGA_LAMBDA_ANALYSIS_BEAM_ENERGY_RGA=10.6 #NOTE: CHANGE AS NEEDED.
export RGA_LAMBDA_ANALYSIS_TPID_RGA=2212 #NOTE: CHANGE AS NEEDED.

# Set prexisting HIPO data paths for CLAS12-Analysis #NOTE: CHANGE AS NEEDED
export MC_RGA_DIR="/cache/clas12/rg-a/production/montecarlo/clasdis_pass1/fall2018/torus+1/v1/bkg50nA_10604MeV/*.hipo"
export DT_RGA_DIR="/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/nSidis/"

# Set numeric values for saga
export RGA_LAMBDA_ANALYSIS_PDG_ALPHA=0.747 #NOTE: This is the lambda decay asymmetry parameter from the PDG
export DT_RGA_BEAM_POLARIZATION=0.8922 #NOTE: This is the luminosity averaged beam polarization
export MC_RGA_BEAM_POLARIZATION=1.0

# Set variables for clas12-config location
export RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME="/work/clas12/users/$USER/clas12-config"

# Set output directory for slurm job stdout and stderr
export RGA_LAMBDA_ANALYSIS_FARM_OUT="/farm_out/$USER"

# Set directory paths
export RGA_LAMBDA_ANALYSIS_VOL_DIR="/volatile/clas12/users/$USER/rga_lambda_analysis" #NOTE: CHANGE AS NEEDED.

# Set image paths
export RGA_LAMBDA_ANALYSIS_GEMC_IMG="gemc_dev-almalinux94/"
export RGA_LAMBDA_ANALYSIS_CCFA_IMG="analysis_latest.sif"
export RGA_LAMBDA_ANALYSIS_C12A_IMG="clas12-analysis.sif"
export RGA_LAMBDA_ANALYSIS_SAGA_IMG="saga.sif"

# Set the HPC partition on which you wish to run
export RGA_LAMBDA_ANALYSIS_HPC_PARTITION="production"
export RGA_LAMBDA_ANALYSIS_HPC_ACCOUNT="clas12"

# Set gemc variables
export RGA_LAMBDA_ANALYSIS_GEMC_VERSION="5.10"
export RGA_LAMBDA_ANALYSIS_COATJAVA_VERSION="11.1.0"

#----- LOAD VARIABLES -----#

# Load and overwrite variables from env.txt
if [ -f env.txt ]; then
    # ignore lines starting with # and blank lines
    export $(grep -v '^#' env.txt | xargs)
fi

#----- STATIC VARIABLES -----#

# Set variables for this project
export RGA_LAMBDA_ANALYSIS_HOME="$PWD"

#----- DEPENDENT VARIABLES -----#

# Set gcard and yaml files for simulation
export RGA_LAMBDA_ANALYSIS_GCARD_RGA_GEMC="$RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME/gemc/$RGA_LAMBDA_ANALYSIS_GEMC_VERSION/rga_fall2018.gcard"
export RGA_LAMBDA_ANALYSIS_YAML_RGA_COAT="$RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME/coatjava/$RGA_LAMBDA_ANALYSIS_COATJAVA_VERSION/rga_fall2018.yaml"
export RGA_LAMBDA_ANALYSIS_ROOT_FILES_DIR="$RGA_LAMBDA_ANALYSIS_VOL_DIR/root_files/mc_rga_sss"
export RGA_LAMBDA_ANALYSIS_NFILES=$(( ($RGA_LAMBDA_ANALYSIS_NEVENTS + $RGA_LAMBDA_ANALYSIS_NMAX - 1) / $RGA_LAMBDA_ANALYSIS_NMAX )) # Calculate number of files needed

# Set project HIPO data paths for CLAS12-Analysis
export MC_RGA_SSS_DIR="$RGA_LAMBDA_ANALYSIS_VOL_DIR/jobs/simulation/mc_rga_sss/dst"

# Set command for gemc
RGA_LAMBDA_ANALYSIS_GEMC_COMMAND() {
    apptainer exec -B $RGA_LAMBDA_ANALYSIS_VOL_DIR,$RGA_LAMBDA_ANALYSIS_HOME,$RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME $RGA_LAMBDA_ANALYSIS_GEMC_IMG \
    bash -c "module use /cvmfs/oasis.opensciencegrid.org/jlab/geant4/modules; \
    module load gemc/$RGA_LAMBDA_ANALYSIS_GEMC_VERSION; \
    export GEMC_DATA_DIR=/cvmfs/oasis.opensciencegrid.org/jlab/geant4/almalinux9-gcc11/clas12Tags/$RGA_LAMBDA_ANALYSIS_GEMC_VERSION; \
    gemc $@"
}
export -f RGA_LAMBDA_ANALYSIS_GEMC_COMMAND

# Set variables for clas12 container forge analysis image
RGA_LAMBDA_ANALYSIS_RECON_UTIL_COMMAND() {
    apptainer exec -B $RGA_LAMBDA_ANALYSIS_VOL_DIR,$RGA_LAMBDA_ANALYSIS_HOME,$RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME $RGA_LAMBDA_ANALYSIS_CCFA_IMG \
    bash -c "/opt/coatjava/bin/recon-util $@"
}
export -f RGA_LAMBDA_ANALYSIS_RECON_UTIL_COMMAND
export RGA_LAMBDA_ANALYSIS_HIPO_UTILS_COMMAND="apptainer exec -B \
$DT_RGA_DIR,$MC_RGA_DIR,$MC_RGA_SSS_DIR,$RGA_LAMBDA_ANALYSIS_VOL_DIR,$RGA_LAMBDA_ANALYSIS_HOME,$RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME $RGA_LAMBDA_ANALYSIS_CCFA_IMG \
bash /opt/coatjava/bin/hipo-utils"

# Set variables for clas12-analysis
export RGA_LAMBDA_ANALYSIS_C12ANALYSIS_COMMAND="apptainer run -B \
$DT_RGA_DIR,$MC_RGA_DIR,$MC_RGA_SSS_DIR,$RGA_LAMBDA_ANALYSIS_VOL_DIR,$RGA_LAMBDA_ANALYSIS_HOME $RGA_LAMBDA_ANALYSIS_C12A_IMG"

# Set variables for saga
RGA_LAMBDA_ANALYSIS_SAGA_COMMAND() {
    apptainer exec -B \
    $RGA_LAMBDA_ANALYSIS_VOL_DIR,$RGA_LAMBDA_ANALYSIS_HOME $RGA_LAMBDA_ANALYSIS_SAGA_IMG \
    bash -c "$@"
}
export -f RGA_LAMBDA_ANALYSIS_SAGA_COMMAND
