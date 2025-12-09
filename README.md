# RGA $\Lambda$ Analysis

This is a storage repository for macros and scripts used run the RGA $\Lambda$ hyperon
spin transfer analysis.  It also allows one to simulate and process [string spinner] generated events in CLAS12.

# Prerequisites
* Assumedly, you are working on ifarm and can use slurm to submit jobs
* [clas12-config](https://github.com/JeffersonLab/clas12-config.git) (`git clone` this and set the path manually in `env.txt`)
* **Make sure you set the torus field to -1 in the gcard!**
* [`gemc`](https://github.com/gemc)
* [`clas12 container forge analysis`](https://pages.jlab.org/hallb/clas12/container-forge/)
* [`clas12-Analysis`](https://github.com/mfmceneaney/CLAS12-Analysis.git)
* [`saga`](https://github.com/mfmceneaney/saga.git)

From the prerequisites above you should check for system installations of the following on ifarm:
* `root`
* `gemc`
* `recon-util` (from clas12 container forge)
* `hipo-utils` (from clas12 container forge)

Otherwise install images with either singularity or apptainer.  Note that you may need to set
the cache and tmp directories for these to some directory capable of housing large files.
For example, on the Duke Compute cluster add the following to your startup script.
```bash
# Set container cache and tmp directory to cwork
export CWORK_DIR=/cwork/$USER/
export APPTAINER_CACHEDIR=$CWORK_DIR
export APPTAINER_TMPDIR=$CWORK_DIR
export SINGULARITY_CACHEDIR=$CWORK_DIR
export SINGULARITY_TMPDIR=$CWORK_DIR
```

Here we will use apptainer to install the necessary images.  You will need to set the path to each image in `env.txt`.

Install `gemc`:
```bash
apptainer pull gemc_dev-almalinux94/ docker://jeffersonlab/gemc:dev-almalinux94
```

Install `clas12 container forge analysis`:
```bash
apptainer pull docker://codecr.jlab.org/hallb/clas12/container-forge/analysis:latest
```

Install `clas12-analysis`:
```bash
apptainer pull clas12-analysis.sif oras://ghcr.io/mfmceneaney/clas12-analysis:latest
```

Install `saga` (or another container that will run CERN `root`):
```bash
apptainer pull saga.sif oras://ghcr.io/mfmceneaney/saga:latest
```

For running the python scripts in [pyscripts](pyscripts),
you will need the saga python modules.  Create a python virtual environment and install the base dependencies.
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```
Then clone the `saga` repository and install saga.
```bash
git clone https://github.com/mfmceneaney/saga.git
cd saga
pip install -e .
```

# Installation

Begin by cloning the repository:
```bash
git clone https://github.com/mfmceneaney/rga_lambda_analysis.git
```

Update the paths and commands used in the environment script&mdash;[bin/env.sh](bin/env.sh) or [bin/env.csh](bin/env.csh)&mdash;by creating a file `env.txt` in the root of this repository.
In this file you will need to manually set variables used in the environment script depending on your local installation paths and the paths for existing data and MC samples you wish to use:
`RGA_LAMBDA_ANALYSIS_VOL_DIR`, `RGA_LAMBDA_ANALYSIS_HOME`, `RGA_LAMBDA_ANALYSIS_*_IMG`, etc.

After configuring your environment file, source the environment and run the setup script.
```bash
source bin/env.sh
./bin/setup.sh
```

Then add the following to your (bash) startup script:
```bash
# Set up rga_lambda_analysis projections https://github.com/mfmceneaney/rga_lambda_analysis.git
pushd /path/to/rga_lambda_analysis >> /dev/null
source bin/env.sh
popd >> /dev/null
```

To run the python download scripts you will also need to create a virtual environment and install the dependencies.
```bash
python3 -m venv /path/to/new/venv
source /path/to/new/venv/bin/activate
pip install -r $RGA_LAMBDA_ANALYSIS_HOME/requirements.txt
/path/to/new/venv/bin/playwright install
```

# Overview

## String Spinner Simulation

First, you must produce simulation HIPO files.
You will need a set of ROOT file containing a tree of ordered, flattened events.
Each entry in the tree should include the following entries:
- `iEvent` : event index \[0,inf\]
- `iPos` : event position index : \[1,inf\]
- `status` : Lund particle status
- `id` : Lund particle PID
- `iMother1` : mother particle Lund index : \[1,inf\]
- `iDaughter1` : first daughter particle Lund index : \[1,inf\]
- `Px` : x-momentum \[GeV\]
- `Py` : y-momentum \[GeV\]
- `Pz` : z-momentum \[GeV\]
- `E` : energy \[GeV\]

Note that the mass will be set automatically using energy conservation, and the charge will also be set automatically based on the PID.

:red_circle: The first three entries (`iPos`$\in[1,2,3]$) of each event should be the incoming electron, the target proton, and the outgoing electron.  The virtual photon will be inferred from the incoming and outgoing electrons using momentum and energy conservation and inserted at index `iPos`$=3$.

The Lund particle status of the first three particles as well as any particle with a non-zero daughter index will be set to `21` so that they are not considered by gemc, while the remaining particles, including the scattered electron, will be assumed to all be final state particles.  Hence, their Lund particle status will be set to `1`.

To download individual ROOT files containing the String Spinner events, run:
```bash
$RGA_LAMBDA_ANALYSIS_HOME/bin/download_root_file.sh <DONWLOAD_URL> <FILENAME>
```

To download a set of ROOT files from a regular expression and cernbox url use the python script instead, assuming your python virtual environment is active.  To see the usage, run:
```bash
python3 $RGA_LAMBDA_ANALYSIS_HOME/pyscripts/cernbox_downloader.py --help
```

Put all the downloaded files in a convenient directory, e.g., `$RGA_LAMBDA_ANALYSIS_VOL_DIR/root_files/mc_rga`, and make `$RGA_LAMBDA_ANALYSIS_ROOT_FILES_DIR` point to this in the environment script.

Then, run the simulation jobs.
```bash
cd $RGA_LAMBDA_ANALYSIS_HOME/jobs/simulation/mc_rga
touch jobs.txt
./setup.sh >> jobs.txt
cd -
```
These may take a while.

## Event Selection

Then, you have to produce channel-specific event-level ROOT files using the directories in `jobs/c12analysis/`.
Make sure to update the paths to existing simulation and data directories for, e.g. MC RGA or Data RGA, in your environment script.
To submit these jobs, cd into the relevant directory and run:
```bash
touch jobs.txt
./setup.sh >> jobs.txt
```

Once these jobs finish, the root files will now be accessible in the directories:
- `$RGA_LAMBDA_ANALYSIS_VOL_DIR/jobs/c12analysis/<run group>`

## Asymmetry Extraction

Run kinematics jobs by going into each directory and manually submitting:
```bash
for file in jobs/saga/test_getBinKinematics*; do
    echo $file
    cd $file
    touch jobs.txt
    ./setup.sh >> jobs.txt
    cd -
    echo
done
```
Then, run the asymmetry injection and extraction studies using the `pyscripts/orchestrate*.py` files.

Finally, aggregate results from injection studies, compute systematics, and plot kinematics and bin schemes with the remaining scripts in `pyscripts`.

#

Contact: matthew.mceneaney@duke.edu
