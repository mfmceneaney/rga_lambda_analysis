# Basic imports
import os
import numpy as np

# Import saga modules
from saga.orchestrate import create_jobs, submit_jobs

# Set dry run to `False` once you are sure you want to submit.
dry_run=True

# Set base directories
run_groups = ["mc_rga","mc_rga_sss"]
channels = ["ppim"]
base_dirs = [
    os.path.abspath(
        os.path.join(
            os.environ['RGA_LAMBDA_ANALYSIS'],
            "jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/"
        )
    ) for rg in run_groups for ch in channels
]

# Loop base directories for basic asymmetry injections
for base_dir in base_dirs:

    # Create job submission structure
    asyms = [-0.1,-0.01,0.0,0.01,0.1]
    asymfitvars = {"asymfitvars":["costheta1","costheta2","costhetaT","costhetaTy"]}
    sgasyms = {"sgasyms":[[a1] for a1 in asyms]}
    bgasyms = {"bgasyms":[[a1] for a1 in asyms]}
    seeds   = {"inject_seed":[2**i for i in range(16)]}

    # Set job file paths and configs
    submit_path =  os.path.join(base_dir,"submit.sh")
    yaml_path   =  os.path.join(base_dir,"args.yaml")
    out_path    =  os.path.join(base_dir,"jobs.txt")
    configs = dict(
        sgasyms,
        **bgasyms,
        **asymfitvars,
        **seeds
    )

    # Create job directories and submit jobs
    create_jobs(configs,base_dir,submit_path,yaml_path)
    submit_jobs(configs,base_dir,submit_path,out_path,dry_run=dry_run)

# Loop base directories and inject extra signal asymmetry
for base_dir in base_dirs:

    # Create job submission structure
    asyms = [-0.1,-0.01,0.0,0.01,0.1]
    sgasyms = {"sgasyms":[[a1,a2] for a1 in asyms for a2 in asyms]}
    asymfitvars = ["costheta1","costheta2","costhetaT","costhetaTy"]
    asymfitvars = {"asymfitvars":[[afv,"phi_h_ppim"] for afv in asymfitvars]}
    seeds   = {"inject_seed":[2**i for i in range(16)]}
    fsgasyms_xs_pu_formula = f"(float)({os.environ['RGA_LAMBDA_ANALYSIS_PDG_ALPHA']}*depol1*sgasyms[0]*asymfitvars[0])"
    fsgasyms = {
        "fsgasyms_xs_pu_formula": [fsgasyms_xs_pu_formula],
    }
    aliases  = {
        "fsgasyms_xs_pu_formula":{
            str(el):"2ASYM"
            for el in fsgasyms["fsgasyms_xs_pu_formula"]
        },
    }

    # Set replacements for default yaml keys
    asymfitvar_lims = {"asymfitvar_titles": ["cos(#theta_{LL'})", "#phi_{p#pi^{-}}"]}
    asymfitvar_lims = {"asymfitvar_lims": [[-1.0, 1.0], [0.0, 6.28]]}
    asymfitpar_inits = {"asymfitpar_inits": [0.0, 0.0]}
    asymfitpar_inits = {"asymfitpar_lims": [[-0.5, 0.5], [-0.5, 0.5]]}
    replacements = dict(
        asymfitvar_lims,
        **asymfitpar_inits,
        **asymfitpar_inits,
    )

    # Set job file paths and configs
    submit_path =  os.path.join(base_dir,"submit.sh")
    yaml_path   =  os.path.join(base_dir,"args.yaml")
    out_path    =  os.path.join(base_dir,"jobs.txt")
    configs = dict(
        sgasyms,
        **bgasyms,
        **seeds
    )

    # Create job directories and submit jobs
    create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases,replacements=replacements)
    submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)
