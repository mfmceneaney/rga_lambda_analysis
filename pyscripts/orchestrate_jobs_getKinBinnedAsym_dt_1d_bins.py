# Basic imports
import os

# Import saga modules
from saga.orchestrate import create_jobs, submit_jobs
from saga.data import load_yaml

# Set dry run to `False` once you are sure you want to submit.
dry_run=True

# Set base directories
run_groups = ["dt_rga"]
channels = ["ppim"]
base_dirs = [
    os.path.abspath(
        os.path.join(
            os.environ['RGA_LAMBDA_ANALYSIS'],
            "jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/"
        )
    ) for rg in run_groups for ch in channels
]

# Set paths for 1D bin scheme yaml for splitting
YAML_DIR = os.path.abspath(os.path.join(os.environ['RGA_LAMBDA_ANALYSIS'],'yamls'))
yaml_paths = [
    os.path.join(YAML_DIR,f'out_1d_bins_{ch}.yaml') for rg in run_groups for ch in channels
]

# Loop base directories with crystal ball signal
for base_dir, yaml_path in zip(base_dirs,yaml_paths):

    # Create job submission structure

    # Split binschemes with aliases
    binschemes = load_yaml(yaml_path)
    binschemes = {"binschemes":[{el:binschemes[el]} for el in binschemes]}
    aliases    = {"binschemes":{
                        str(el):list(el.keys())[0]+"_binscheme"
                        for el in binschemes["binschemes"]
                    }
                }

    # Set job file paths and configs
    submit_path = os.path.join(base_dir,"submit.sh")
    yaml_path   = os.path.join(base_dir,"args.yaml")
    out_path    = os.path.join(base_dir,"jobs.txt")
    configs = dict(
        binschemes,
    )

    # Create job directories and submit jobs
    create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases)
    submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)

# Create dictionary of mass fit yaml file path dictionary
massfit_type = "dt_massfit_gs"
massfit_yamlfile_maps = { "massfit_yamlfile_map":{
        "scheme_Q2_bin_0": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_Q2_bin_0.yaml",
        "scheme_Q2_bin_1": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_Q2_bin_1.yaml",
        "scheme_Q2_bin_2": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_Q2_bin_2.yaml",
        "scheme_Q2_bin_3": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_Q2_bin_3.yaml",
        "scheme_Q2_bin_4": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_Q2_bin_4.yaml",
        "scheme_W_bin_0": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_W_bin_0.yaml",
        "scheme_W_bin_1": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_W_bin_1.yaml",
        "scheme_W_bin_2": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_W_bin_2.yaml",
        "scheme_W_bin_3": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_W_bin_3.yaml",
        "scheme_W_bin_4": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_W_bin_4.yaml",
        "scheme_y_bin_0": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_y_bin_0.yaml",
        "scheme_y_bin_1": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_y_bin_1.yaml",
        "scheme_y_bin_2": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_y_bin_2.yaml",
        "scheme_y_bin_3": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_y_bin_3.yaml",
        "scheme_y_bin_4": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_y_bin_4.yaml",
        "scheme_x_bin_0": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_x_bin_0.yaml",
        "scheme_x_bin_1": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_x_bin_1.yaml",
        "scheme_x_bin_2": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_x_bin_2.yaml",
        "scheme_x_bin_3": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_x_bin_3.yaml",
        "scheme_x_bin_4": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_x_bin_4.yaml",
        "scheme_xF_ppim_bin_0": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_xF_ppim_bin_0.yaml",
        "scheme_xF_ppim_bin_1": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_xF_ppim_bin_1.yaml",
        "scheme_xF_ppim_bin_2": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_xF_ppim_bin_2.yaml",
        "scheme_xF_ppim_bin_3": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_xF_ppim_bin_3.yaml",
        "scheme_xF_ppim_bin_4": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_xF_ppim_bin_4.yaml",
        "scheme_z_ppim_bin_0": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_z_ppim_bin_0.yaml",
        "scheme_z_ppim_bin_1": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_z_ppim_bin_1.yaml",
        "scheme_z_ppim_bin_2": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_z_ppim_bin_2.yaml",
        "scheme_z_ppim_bin_3": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_z_ppim_bin_3.yaml",
        "scheme_z_ppim_bin_4": f"/{RGA_LAMBDA_ANALYSIS_HOME}/yamls/{massfit_type}/scheme_z_ppim_bin_4.yaml",
    }
}

# Loop base directories with gaussian signal
for base_dir, yaml_path in zip(base_dirs,yaml_paths):

    # Create job submission structure

    # Split binschemes with aliases
    binschemes = load_yaml(yaml_path)
    binschemes = {"binschemes":[{el:binschemes[el]} for el in binschemes]}
    massfit_yamlfile_maps = {
        "massfit_yamlfile_map": massfit_yamlfile_map
    }
    aliases    = {
        "binschemes":{
            str(el):list(el.keys())[0]+"_binscheme"
            for el in binschemes["binschemes"]
        },
        "massfit_yamlfile_map":{
            "massfit_yamlfile_map":"gs",
        }
    }

    # Set job file paths and configs
    submit_path = os.path.join(base_dir,"submit.sh")
    yaml_path   = os.path.join(base_dir,"args.yaml")
    out_path    = os.path.join(base_dir,"jobs.txt")
    configs = dict(
        binschemes,
        **massfit_yamlfile_maps,
    )

    # Create job directories and submit jobs
    create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases)
    submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)
