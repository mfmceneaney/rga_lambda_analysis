# Basic imports
import os
import argparse

# Import saga modules
from saga.orchestrate import create_jobs, submit_jobs
from saga.data import load_yaml
from saga.aggregate import get_scheme_vars, get_binscheme_cuts_and_ids

# Parse arguments
parser = argparse.ArgumentParser(description='Script to submit `getKinBinnedAsym` and `getKinBinnedHB` jobs on RGA data for the `Lambda -> proton pion` channel')
parser.add_argument('--dry_run', default=True, help='Dry run without job submission')
parser.add_argument('--basic', default=True, help='Submit sideband subtraction asymmetry extraction jobs')
parser.add_argument('--massfit_types', default=None, help='Submit mass fit signal type jobs', nargs="*", choices=["gauss","doublegauss","landau","breitwigner","crystalball"])
parser.add_argument('--cos_phi', default=False, help='Submit cos_phi difference jobs')
parser.add_argument('--no_sbs', default=False, help='Submit no sb subtraction jobs')
parser.add_argument('--splot', default=False, help='Submit sPlot asymmetry extraction jobs')
parser.add_argument('--asymfitvars', default=["costheta1","costheta2","costhetaT","costhetaTy"], help='Lambda decay angle fit variables to use', nargs="+", choices=["costheta1","costheta2","costhetaT","costhetaTy"])
parser.add_argument('--rgs', default=["dt_rga"], help='Run group', nargs="+", choices=["dt_rga"])
args = parser.parse_args()

# Set dry run to `False` once you are sure you want to submit.
dry_run=args.dry_run

# Set base directories
run_groups = args.rgs
channels = ["ppim"]
RGA_LAMBDA_ANALYSIS_HOME = os.environ['RGA_LAMBDA_ANALYSIS_HOME']
YAML_DIR = os.path.abspath(os.path.join(RGA_LAMBDA_ANALYSIS_HOME,'yamls'))

# Loop run groups and channels with crystal ball signal
if args.basic:
    for rg in run_groups:
        for ch in channels:

            # Set base directory and bin scheme yaml
            base_dir = os.path.abspath(
                os.path.join(
                    RGA_LAMBDA_ANALYSIS_HOME,
                    "jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/"
                )
            )
            binscheme_yaml_path = os.path.join(YAML_DIR,f'out_1d_bins_{ch}.yaml')

            # Create job submission structure
            binschemes  = load_yaml(binscheme_yaml_path)
            asymfitvars = {"asymfitvars":args.asymfitvars}
            binschemes  = {"binschemes":[{el:binschemes[el]} for el in binschemes]}
            aliases     = {
                "binschemes":{
                    str(el):list(el.keys())[0]+"_binscheme"
                    for el in binschemes["binschemes"]
                }
            }

            # Set replacements
            replacements = None

            # Set job file paths and configs
            submit_path = os.path.join(base_dir,"submit.sh")
            yaml_path   = os.path.join(base_dir,"args.yaml")
            out_path    = os.path.join(base_dir,"jobs.txt")
            configs = dict(
                asymfitvars,
                **binschemes,
            )

            # Create job directories and submit jobs
            create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases,replacements=replacements)
            submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)

# Loop run groups and channels with gaussian signal
if args.massfit_types is not None:
    for massfit_type in massfit_types:

        # Loop run groups and channels
        for rg in run_groups:
            for ch in channels:

                # Set base directory and bin scheme yaml
                base_dir = os.path.abspath(
                    os.path.join(
                        RGA_LAMBDA_ANALYSIS_HOME,
                        "jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/"
                    )
                )
                binscheme_yaml_path = os.path.join(YAML_DIR,f'out_1d_bins_{ch}.yaml')

                # Create job submission structure
                binschemes  = load_yaml(binscheme_yaml_path)
                asymfitvars = {"asymfitvars":args.asymfitvars}
                binschemes  = {"binschemes":[{el:binschemes[el]} for el in binschemes]}

                # Create dictionary of mass fit yaml file path dictionary
                massfit_yamlfile_map = {
                    f"scheme_{binscheme}_bin_{binid}": \
                    os.path.join(
                        f"{RGA_LAMBDA_ANALYSIS_HOME}/yamls/massfit/{rg}/{massfit_type}/"
                        f"scheme_{binvar}_bin_{binid}.yaml",
                    ) for binid in range(len(get_binscheme_cuts_and_ids(binscheme)[2]))
                    for binscheme in binschemes 
                }
                print("INFO: massfit_yamlfile_map = {")
                for key in massfit_yamlfile_map:
                    print(f"INFO: \t{key}: {massfit_yamlfile_map[key]},")
                print("INFO: }")
                massfit_yamlfile_maps = {
                    "massfit_yamlfile_map": massfit_yamlfile_map
                }

                # Set aliases
                aliases     = {
                    "binschemes":{
                        str(el):list(el.keys())[0]+"_binscheme"
                        for el in binschemes["binschemes"]
                    },
                    "massfit_yamlfile_map":{
                        str(massfit_yamlfile_map):"gs",
                    }
                }

                # Set replacements
                replacements = None

                # Set job file paths and configs
                submit_path = os.path.join(base_dir,"submit.sh")
                yaml_path   = os.path.join(base_dir,"args.yaml")
                out_path    = os.path.join(base_dir,"jobs.txt")
                configs = dict(
                    asymfitvars,
                    **binschemes,
                    **massfit_yamlfile_maps,
                )

                # Create job directories and submit jobs
                create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases,replacements=replacements)
                submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)

# Loop run groups and channels with 2 cos phi regions
if args.cos_phi:
    for rg in run_groups:
        for ch in channels:

            # Set base directory and bin scheme yaml
            base_dir = os.path.abspath(
                os.path.join(
                    RGA_LAMBDA_ANALYSIS_HOME,
                    "jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/"
                )
            )
            binscheme_yaml_path = os.path.join(YAML_DIR,f'out_1d_bins_{ch}.yaml')

            # Create job submission structure with fit variables and cos phi cuts
            asymfitvars = {"asymfitvars":args.asymfitvars}
            args_yaml_path = os.path.join(base_dir,"args.yaml")
            args_yaml = load_yaml(args_yaml_path)
            cuts = args_yaml["cuts"]
            cuts_pos_cos_phi = cuts + " && !(phi_h_ppim<TMath::Pi()/2 || phi_h_ppim>=3*TMath::Pi()/2)"
            cuts_neg_cos_phi = cuts + " && (phi_h_ppim<TMath::Pi()/2 || phi_h_ppim>=3*TMath::Pi()/2)"
            cutss = {"cuts",[cuts_pos_cos_phi,cuts_neg_cos_phi]}
            aliases     = {
                "cuts":{
                    cuts_pos_cos_phi:"pos_cos_phi",
                    cuts_neg_cos_phi:"neg_cos_phi",
                }
            }

            # Set replacements
            replacements = {
                "binschemes_paths": [f"{RGA_LAMBDA_ANALYSIS_HOME}/yamls/out_full_bin_ppim.yaml"],
            }

            # Set job file paths and configs
            submit_path = os.path.join(base_dir,"submit.sh")
            yaml_path   = os.path.join(base_dir,"args.yaml")
            out_path    = os.path.join(base_dir,"jobs.txt")
            configs = dict(
                asymfitvars,
                **cutss,
            )

            # Create job directories and submit jobs
            create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases,replacements=replacements)
            submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)

# Loop run groups and channels with many bins and no sb subtraction binned in mass_ppim
if args.no_sbs:
    for rg in run_groups:
        for ch in channels:

            # Set base directory and bin scheme yaml
            base_dir = os.path.abspath(
                os.path.join(
                    RGA_LAMBDA_ANALYSIS_HOME,
                    "jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/"
                )
            )
            binscheme_yaml_path = os.path.join(YAML_DIR,f'out_1d_bins_{ch}.yaml')

            # Create job submission structure with fit variables and mass_ppim binscheme and no sideband subtraction
            asymfitvars = {"asymfitvars":args.asymfitvars}
            use_sb_subtractions = {"use_sb_subtraction": [False]}
            binscheme_yaml_path = f"{RGA_LAMBDA_ANALYSIS_HOME}/yamls/out_1d_bins_ppim_no_bg_correction.yaml"
            binscheme_yaml = load_yaml(binscheme_yaml_path)["mass_ppim"]
            binschemes = {"binschemes":[binscheme_yaml]}
            aliases     = {
                "binschemes":{
                    str(binscheme_yaml): "mass_ppim_binscheme",
                },
            }

            # Set replacements
            replacements = None

            # Set job file paths and configs
            submit_path = os.path.join(base_dir,"submit.sh")
            yaml_path   = os.path.join(base_dir,"args.yaml")
            out_path    = os.path.join(base_dir,"jobs.txt")
            configs = dict(
                asymfitvars,
                **binschemes,
                **use_sb_subtractions,
            )

            # Create job directories and submit jobs
            create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases,replacements=replacements)
            submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)

# Loop run groups and channels with many bins and no sb subtraction binned in xF_ppim with no xF>0 cut
if args.no_sbs:
    for rg in run_groups:
        for ch in channels:

            # Set base directory and bin scheme yaml
            base_dir = os.path.abspath(
                os.path.join(
                    RGA_LAMBDA_ANALYSIS_HOME,
                    "jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/"
                )
            )
            binscheme_yaml_path = os.path.join(YAML_DIR,f'out_1d_bins_{ch}_no_bg_correction.yaml')

            # Create job submission structure with fit variables and xF_ppim binscheme and no sideband subtraction
            asymfitvars = {"asymfitvars":args.asymfitvars}
            use_sb_subtractions = {"use_sb_subtraction": [False]}
            binscheme_yaml = load_yaml(binscheme_yaml_path)["xF_ppim"]
            binschemes = {"binschemes":[binscheme_yaml]}
            aliases     = {
                "binschemes":{
                    str(binscheme_yaml): "xF_ppim_binscheme",
                },
            }

            # Set replacements
            args_yaml_path = os.path.join(base_dir,"args.yaml")
            args_yaml = load_yaml(args_yaml_path)
            cuts = args_yaml["cuts"].replace(" && xF_ppim>0.0","")
            replacements = {
                "cuts": cuts,
            }

            # Set job file paths and configs
            submit_path = os.path.join(base_dir,"submit.sh")
            yaml_path   = os.path.join(base_dir,"args.yaml")
            out_path    = os.path.join(base_dir,"jobs.txt")
            configs = dict(
                asymfitvars,
                **binschemes,
                **use_sb_subtractions,
            )

            # Create job directories and submit jobs
            create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases,replacements=replacements)
            submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)

# Loop run groups and channels with 1 bin and splot
if args.splot:
    for rg in run_groups:
        for ch in channels:

            # Set base directory and bin scheme yaml
            base_dir = os.path.abspath(
                os.path.join(
                    RGA_LAMBDA_ANALYSIS_HOME,
                    "jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/"
                )
            )
            binscheme_yaml_path = os.path.join(YAML_DIR,f'out_1d_bins_{ch}.yaml')

            # Create job submission structure with fit variables and cos phi cuts
            binschemes  = load_yaml(yaml_path)
            asymfitvars = {"asymfitvars":args.asymfitvars}
            binschemes  = {"binschemes":[{el:binschemes[el]} for el in binschemes]}
            use_splot   = {"use_splot": [True]}
            aliases     = {
                "binschemes":{
                    str(el):list(el.keys())[0]+"_binscheme"
                    for el in binschemes["binschemes"]
                }
            }

            # Set replacements
            replacements = {
                "use_sb_subtraction": False,
            }

            # Set job file paths and configs
            submit_path = os.path.join(base_dir,"submit.sh")
            yaml_path   = os.path.join(base_dir,"args.yaml")
            out_path    = os.path.join(base_dir,"jobs.txt")
            configs = dict(
                asymfitvars,
                **binschemes,
                **use_splot,
            )

            # Create job directories and submit jobs
            create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases,replacements=replacements)
            submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)
