# CLAS12 RGA $\Lambda$ Baryon Polarization Analysis

This is a storage repository for yamls and scripts used to run the CLAS12 RGA $\Lambda$ baryon polarization analysis.

# Prerequisites
* Assumedly, you are working on ifarm and can use slurm to submit jobs
* [CLAS12-Analysis](https://github.com/mfmceneaney/CLAS12-Analysis.git)
* [saga](https://github.com/mfmceneaney/saga.git)

# Installation

Begin by cloning the repository:
```bash
git clone https://github.com/mfmceneaney/rga_lambda_analysis.git
```

Make sure all the paths in the environment script&mdash;[bin/env.sh](bin/env.sh) or [bin/env.csh](bin/env.csh)&mdash;are correct for you.
In particular, you will need to manually set these variables in the environment script depending on your local installation paths and the paths for existing data and MC samples you wish to use:
`RGA_LAMBDA_ANALYSIS_HOME`, `RGA_LAMBDA_ANALYSIS_VOL_DIR`, `SAGA_HOME`, `SAGA_BUILD_DIR`, `RGA_MC_DIR`, `RGA_DT_DIR`.
Yaml paths will be set based on the paths given in the environment script.

After configuring your environment script, add the following to your (bash) startup script:
```bash
# Set up RGA Lambda Analysis https://github.com/mfmceneaney/rga_lambda_analysis.git
pushd /path/to/rga_lambda_analysis >> /dev/null
source bin/env.sh
popd >> /dev/null
```

# Overview

First, you have to produce channel-specific event-level ROOT files using the directories in `jobs/c12analysis/`.
Make sure to update the paths to existing simulation and data directories in your environment script.
To submit these jobs, cd into the relevant directory and run:
```bash
touch jobs.txt
./setup.sh >> jobs.txt
```

Configure yamls for jobs running [saga](https://github.com/mfmceneaney/saga.git) by running `$RGA_LAMBDA_ANALYSIS_HOME/bin/setup.sh`.

Run kinematics (and bin migration) jobs by going into each directory and manually submitting:
```bash
for file in jobs/saga/test_getBinKinematics*; do
    echo $file
    cd $file
    touch jobs.txt
    sbatch $PWD/submit.sh >> jobs.txt
    cd -
    echo
done
```

Run injection studies and data studies using the `pyscripts/orchestrate*.py` files.

Finally, aggregate results from injection and data studies using the `pyscripts/aggregate*.py` files.

#

Contact: matthew.mceneaney@duke.edu
