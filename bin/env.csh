#!/bin/tcsh

#----- DEFAULT VARIABLES -----#

# Set variables for string spinner simulation
setenv RGA_LAMBDA_ANALYSIS_TREE "tree" #NOTE: CHANGE AS NEEDED.
setenv RGA_LAMBDA_ANALYSIS_NEVENTS 10000000 # Number of events in each file #NOTE: CHANGE AS NEEDED.
setenv RGA_LAMBDA_ANALYSIS_NMAX 10000 # Maximum number of events allowed in each output lund file

# Set beam energies and target lund pids for CLAS12-Analysis #NOTE: CHANGE AS NEEDED
setenv RGA_LAMBDA_ANALYSIS_BEAM_ENERGY_RGA 10.6 #NOTE: CHANGE AS NEEDED.
setenv RGA_LAMBDA_ANALYSIS_TPID_RGA 2212 #NOTE: CHANGE AS NEEDED.

# Set prexisting HIPO data paths for CLAS12-Analysis #NOTE: CHANGE AS NEEDED
setenv MC_RGA_DIR "/cache/clas12/rg-a/production/montecarlo/clasdis_pass1/fall2018/torus+1/v1/bkg50nA_10604MeV/*.hipo"
setenv DT_RGA_DIR "/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/nSidis/"

# Set numeric values for saga
setenv RGA_LAMBDA_ANALYSIS_PDG_ALPHA 0.747 #NOTE: This is the lambda decay asymmetry parameter from the PDG
setenv DT_RGA_BEAM_POLARIZATION 0.8922 #NOTE: This is the luminosity averaged beam polarization
setenv MC_RGA_BEAM_POLARIZATION 1.0

# Set variables for clas12-config location
setenv RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME "/work/clas12/users/$USER/clas12-config"

# Set output directory for slurm job stdout and stderr
setenv RGA_LAMBDA_ANALYSIS_FARM_OUT "/farm_out/$USER"

# Set directory paths
setenv RGA_LAMBDA_ANALYSIS_VOL_DIR "/volatile/clas12/users/$USER/rga_lambda_analysis" #NOTE: CHANGE AS NEEDED.

# Set image paths
setenv RGA_LAMBDA_ANALYSIS_GEMC_IMG "gemc_dev-almalinux94/"
setenv RGA_LAMBDA_ANALYSIS_CCFA_IMG "analysis_latest.sif"
setenv RGA_LAMBDA_ANALYSIS_C12A_IMG "clas12-analysis.sif"
setenv RGA_LAMBDA_ANALYSIS_SAGA_IMG "saga.sif"

# Set the HPC partition on which you wish to run
setenv RGA_LAMBDA_ANALYSIS_HPC_PARTITION "production"
setenv RGA_LAMBDA_ANALYSIS_HPC_ACCOUNT "clas12"

# Set gemc variables
setenv RGA_LAMBDA_ANALYSIS_GEMC_VERSION "5.10"
setenv RGA_LAMBDA_ANALYSIS_COATJAVA_VERSION "11.1.0"

#----- LOAD VARIABLES -----#

# Load and overwrite variables from env.txt
if (-e env.txt) then
    foreach line (`grep -v '^#' env.txt`)
        set var = `echo $line | cut -d= -f1`
        set val = `echo $line | cut -d= -f2-`
        setenv $var "$val"
    end
endif

#----- STATIC VARIABLES -----#

# Set variables for this project
setenv RGA_LAMBDA_ANALYSIS_HOME "$PWD"

#----- DEPENDENT VARIABLES -----#

# Set gcard and yaml files for simulation
setenv RGA_LAMBDA_ANALYSIS_GCARD_RGA_GEMC "$RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME/gemc/$RGA_LAMBDA_ANALYSIS_GEMC_VERSION/rga_fall2018.gcard"
setenv RGA_LAMBDA_ANALYSIS_YAML_RGA_COAT "$RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME/coatjava/$RGA_LAMBDA_ANALYSIS_COATJAVA_VERSION/rga_fall2018.yaml"
setenv RGA_LAMBDA_ANALYSIS_ROOT_FILES_DIR "$RGA_LAMBDA_ANALYSIS_VOL_DIR/root_files/mc_rga_sss"

# Calculate number of files needed: (RGA_LAMBDA_ANALYSIS_NEVENTS + RGA_LAMBDA_ANALYSIS_NMAX - 1) / RGA_LAMBDA_ANALYSIS_NMAX
@ tmp = ( $RGA_LAMBDA_ANALYSIS_NEVENTS + $RGA_LAMBDA_ANALYSIS_NMAX - 1 )
@ RGA_LAMBDA_ANALYSIS_NFILES = $tmp / $RGA_LAMBDA_ANALYSIS_NMAX
setenv RGA_LAMBDA_ANALYSIS_NFILES $RGA_LAMBDA_ANALYSIS_NFILES

# Set project HIPO data paths for CLAS12-Analysis
setenv MC_RGA_SSS_DIR "$RGA_LAMBDA_ANALYSIS_VOL_DIR/jobs/simulation/mc_rga_sss/dst"

# Define helper commands as aliases that call sh -c to preserve original bash behavior
alias RGA_LAMBDA_ANALYSIS_GEMC_COMMAND 'sh -c "apptainer exec -B \$RGA_LAMBDA_ANALYSIS_VOL_DIR,\$RGA_LAMBDA_ANALYSIS_HOME,\$RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME \$RGA_LAMBDA_ANALYSIS_GEMC_IMG bash -c \"module use /cvmfs/oasis.opensciencegrid.org/jlab/geant4/modules; module load gemc/\$RGA_LAMBDA_ANALYSIS_GEMC_VERSION; setenv GEMC_DATA_DIR=/cvmfs/oasis.opensciencegrid.org/jlab/geant4/almalinux9-gcc11/clas12Tags/\$RGA_LAMBDA_ANALYSIS_GEMC_VERSION; gemc \" \!*"'

alias RGA_LAMBDA_ANALYSIS_RECON_UTIL_COMMAND 'sh -c "apptainer exec -B \$RGA_LAMBDA_ANALYSIS_VOL_DIR,\$RGA_LAMBDA_ANALYSIS_HOME,\$RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME \$RGA_LAMBDA_ANALYSIS_CCFA_IMG bash -c \"/opt/coatjava/bin/recon-util \" \!*"'

alias RGA_LAMBDA_ANALYSIS_HIPO_UTILS_COMMAND 'sh -c "apptainer exec -B \$DT_RGA_DIR,\$MC_RGA_DIR,\$MC_RGA_SSS_DIR,\$RGA_LAMBDA_ANALYSIS_VOL_DIR,\$RGA_LAMBDA_ANALYSIS_HOME,\$RGA_LAMBDA_ANALYSIS_C12_CONFIG_HOME \$RGA_LAMBDA_ANALYSIS_CCFA_IMG bash /opt/coatjava/bin/hipo-utils"'

alias RGA_LAMBDA_ANALYSIS_C12ANALYSIS_COMMAND 'sh -c "apptainer run -B \$DT_RGA_DIR,\$MC_RGA_DIR,\$MC_RGA_SSS_DIR,\$RGA_LAMBDA_ANALYSIS_VOL_DIR,\$RGA_LAMBDA_ANALYSIS_HOME \$RGA_LAMBDA_ANALYSIS_C12A_IMG"'

alias RGA_LAMBDA_ANALYSIS_SAGA_COMMAND 'sh -c "apptainer exec -B \$RGA_LAMBDA_ANALYSIS_VOL_DIR,\$RGA_LAMBDA_ANALYSIS_HOME \$RGA_LAMBDA_ANALYSIS_SAGA_IMG bash -c \"\!*\""'

