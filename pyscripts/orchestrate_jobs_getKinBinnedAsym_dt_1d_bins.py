# Basic imports
import sys
import os

# Import saga modules
SAGA_HOME = os.environ['SAGA_HOME']
sys.path.append(os.path.abspath(os.path.join(SAGA_HOME,'py')))
from saga.orchestrate import create_jobs, submit_jobs
import saga.aggregate as sagas

# Set dry run to `False` once you are sure you want to submit.
dry_run=True

# Set base directories
run_groups = ["dt_rga"]
channels = ["ppim"]
base_dirs = [
    os.path.abspath(os.path.join(os.environ['RGA_LAMBDA_ANALYSIS_HOME'],"jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/")) for rg in run_groups for ch in channels
]

# Set paths for 1D bin scheme yaml for splitting
YAML_DIR = os.path.abspath(os.path.join(os.environ['RGA_LAMBDA_ANALYSIS_HOME'],'yamls'))
yaml_paths = [
    os.path.join(YAML_DIR,f'out_1d_bins_{ch}.yaml') for rg in run_groups for ch in channels
]

# Loop base directories
for base_dir, yaml_path in zip(base_dirs,yaml_paths):

    # Create job submission structure
    asyms = [0.0]
    sgasyms = {"sgasyms":[[a1,a2,a3] for a1 in asyms for a2 in asyms for a3 in asyms]}
    seeds   = {"inject_seed":[2**i for i in range(1)]}

    # Set job file paths and configs
    submit_path = os.path.join(base_dir,"submit.sh")
    yaml_path   = os.path.join(base_dir,"args.yaml")
    out_path    = os.path.join(base_dir,"jobs.txt")
    configs = dict(
        sgasyms,
        **seeds
    )

    # Create job directories and submit jobs
    create_jobs(configs,base_dir,submit_path,yaml_path)
    submit_jobs(configs,base_dir,submit_path,out_path,dry_run=dry_run)
