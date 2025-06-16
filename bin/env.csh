#!/bin/tcsh

# Set variables for this project
setenv RGA_LAMBDA_ANALYSIS_HOME $PWD
setenv RGH_LAMBDA_ANALYSIS_VOL_DIR "/volatile/clas12/users/$USER/rga_lambda_analysis" #NOTE: CHANGE AS NEEDED.

# Set variables for saga
setenv SAGA_HOME "/work/clas12/users/$USER/saga" #NOTE: CHANGE AS NEEDED. This is a path to your local repo of: https://github.com/mfmceneaney/saga
setenv SAGA_BUILD_DIR "$SAGA_HOME/build"

# Set beam energies and target lund pids for clasdis and CLAS12-Analysis #NOTE: CHANGE AS NEEDED
setenv BEAM_ENERGY_RGA 10.6
setenv TPID_RGA 2212
setenv BPOL_RGA 0.8922 #NOTE: This is for outbending F18 data only.

# Set prexisting HIPO data paths for CLAS12-Analysis #NOTE: CHANGE AS NEEDED
setenv RGA_MC_DIR "/cache/clas12/rg-a/production/montecarlo/clasdis_pass1/fall2018/torus+1/v1/bkg50nA_10604MeV/o"
setenv RGA_DT_DIR "/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/nSidis"
