# Basic imports
import os
import argparse

# Import saga modules
from saga.orchestrate import create_jobs, submit_jobs

# Parse arguments
parser = argparse.ArgumentParser(description='Script to submit `getKinBinnedAsym` and `getKinBinnedHB` jobs on RGA MC for the `Lambda -> proton pion` channel')
parser.add_argument('--dry_run', default=True, help='Dry run without job submission')
parser.add_argument('--splot', default=False, help='Submit splot asymmetry extraction jobs')
parser.add_argument('--massfit_types', default=None, help='Submit mass fit signal type jobs', nargs="*", choices=["gauss","doublegauss","landau","breitwigner","crystalball"])
parser.add_argument('--cos_phi', default=False, help='Submit cos_phi difference jobs')
parser.add_argument('--n_inject_seeds', default=16, help='Number of random injection seeds to use')
parser.add_argument('--asymfitvars', default=["costheta1","costheta2","costhetaT","costhetaTy"], help='Lambda decay angle fit variables to use', nargs="+", choices=["costheta1","costheta2","costhetaT","costhetaTy"])
parser.add_argument('--sgasyms', default=None, help='Signal asymmetries to inject', nargs="*", type=float)
parser.add_argument('--bgasyms', default=None, help='Background asymmetries to inject', nargs="*", type=float)
parser.add_argument('--sgasyms2', default=None, help='Signal asymmetries to inject for additional cos_phi dependent signal term', nargs="*", type=float)
parser.add_argument('--rgs', default=["mc_rga"], help='Run group', nargs="+", choices=["mc_rga","mc_rga_sss"])
args = parser.parse_args()

# Set dry run to `False` once you are sure you want to submit.
dry_run=args.dry_run

# Set base directories
run_groups = args.rgs
channels = ["ppim"]
RGA_LAMBDA_ANALYSIS_HOME = os.environ['RGA_LAMBDA_ANALYSIS_HOME']
YAML_DIR = os.path.abspath(os.path.join(RGA_LAMBDA_ANALYSIS_HOME,'yamls'))

# Loop run groups and channels for basic asymmetry injections
if sgasyms and bgasyms:
    for rg in run_groups:
        for ch in channels:

            # Grab the base directory
            base_dir = os.path.abspath(
                os.path.join(
                    RGA_LAMBDA_ANALYSIS_HOME,
                    "jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/"
                )
            ) 

            # Create job submission structure
            asymfitvars = {"asymfitvars":args.asymfitvars}
            sgasyms = {"sgasyms":[[a1] for a1 in args.sgasyms]}
            bgasyms = {"bgasyms":[[a1] for a1 in args.bgasyms]}
            seeds   = {"inject_seed":[2**i for i in range(args.n_inject_seeds)]}
            aliases = None

            # Set replacements
            replacements = None

            # Set job file paths and configs
            submit_path =  os.path.join(base_dir,"submit.sh")
            yaml_path   =  os.path.join(base_dir,"args.yaml")
            out_path    =  os.path.join(base_dir,"jobs.txt")
            configs = dict(
                asymfitvars,
                **sgasyms,
                **bgasyms,
                **seeds
            )
            if args.splot:
                splot = {"use_splot":True}
                config.update(splot)

            # Create job directories and submit jobs
            create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases,replacements=replacements)
            submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)

# Loop run groups and channels and inject extra signal asymmetry
if sgasyms and sgasyms2:
    for rg in run_groups:
        for ch in channels:

            # Grab the base directory
            base_dir = os.path.abspath(
                os.path.join(
                    RGA_LAMBDA_ANALYSIS_HOME,
                    "jobs/saga/",f"test_getKinBinnedAsym__{rg}__{ch}__1D/"
                )
            ) 

            # Create job submission structure
            sgasyms = {"sgasyms":[[a1,a2] for a1 in args.sgasyms for a2 in args.sgasyms2]}
            asymfitvars = args.asymfitvars
            asymfitvars = {"asymfitvars":[[afv,"phi_h_ppim"] for afv in asymfitvars]}
            seeds   = {"inject_seed":[2**i for i in range(args.n_inject_seeds)]}
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
                asymfitvars,
                **fsgasyms,
                **sgasyms,
                **seeds
            )
            if args.splot:
                splot = {"use_splot":True}
                config.update(splot)

            # Create job directories and submit jobs
            create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases,replacements=replacements)
            submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)

# Loop mass fits with signal pdf type
if args.massfit_types is not None:
    for massfit_type in massfit_types:

        # Loop run groups and channels and inject extra signal asymmetry
        for rg in run_groups:
            for ch in channels:

                # Grab the base directory
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
                sgasyms = {"sgasyms":[[a1] for a1 in args.sgasyms]}
                bgasyms = {"bgasyms":[[a1] for a1 in args.bgasyms]}
                seeds   = {"inject_seed":[2**i for i in range(args.n_inject_seeds)]}
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
                    **sgasyms,
                    **bgasyms,
                    **seeds,
                    **massfit_yamlfile_maps,
                )
                if args.splot:
                    splot = {"use_splot":True}
                    config.update(splot)

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

            # Create job submission structure with fit variables and cos phi cuts
            asymfitvars = {"asymfitvars":args.asymfitvars}
            sgasyms = {"sgasyms":[[a1] for a1 in args.sgasyms]}
            bgasyms = {"bgasyms":[[a1] for a1 in args.bgasyms]}
            seeds   = {"inject_seed":[2**i for i in range(args.n_inject_seeds)]}
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
                **sgasyms,
                **bgasyms,
                **seeds,
            )
            if args.splot:
                splot = {"use_splot":True}
                config.update(splot)

            # Create job directories and submit jobs
            create_jobs(configs,base_dir,submit_path,yaml_path,aliases=aliases,replacements=replacements)
            submit_jobs(configs,base_dir,submit_path,out_path,aliases=aliases,dry_run=dry_run)
