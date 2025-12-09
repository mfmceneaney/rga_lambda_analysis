# Basic imports
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import sys

# Import saga modules
import saga.aggregate as sagas
from saga.aggregate import get_binscheme_cuts_and_ids
from saga.data import load_yaml, load_csv, save_bin_mig_mat_to_csv
from saga.plot import set_default_plt_settings, plot_results

# Parse arguments
parser = argparse.ArgumentParser(description='Script to aggregate `getKinBinnedAsym` and `getKinBinnedHB` jobs on RGA data for the `Lambda -> proton pion` channel')
parser.add_argument('--basic', action="store_true", help='Aggregate sideband subtraction asymmetry extraction jobs')
parser.add_argument('--massfit_types', default=None, help='Aggregate mass fit signal type jobs', nargs="*", choices=["gaus","doublegaus","landau","breitwigner","crystalball"])
parser.add_argument('--cos_phi', action="store_true", help='Aggregate cos_phi difference jobs')
parser.add_argument('--no_sbs_mass_ppim', action="store_true", help='Aggregate no sb subtraction jobs binned in mass_ppim')
parser.add_argument('--no_sbs_xF_ppim', action="store_true", help='Aggregate no sb subtraction jobs binned in xF_ppim')
parser.add_argument('--splot', action="store_true", help='Aggregate sPlot asymmetry extraction jobs')
parser.add_argument('--asymfitvars', default=["costheta1","costheta2","costhetaT","costhetaTy"], help='Lambda decay angle fit variables to use', nargs="+", choices=["costheta1","costheta2","costhetaT","costhetaTy"])
parser.add_argument('--rgs', default=["dt_rga"], help='Run group', nargs="+", choices=["dt_rga"])
parser.add_argument('--methods', default=["HB"], help='Asymmetry extraction method', nargs="+", choices=["HB","Asym"])
parser.add_argument('--use_bin_mig', action="store_true", help='Use bin migration correction')
args = parser.parse_args()

# Set base directory from environment
RGA_LAMBDA_ANALYSIS_HOME = os.environ['RGA_LAMBDA_ANALYSIS_HOME']

# Set up chaining for batched data (specifically `old_dat_path`)
nbatch = 1
nbatches = {"nbatches":[nbatch]}
ibatches = {"ibatch":[i for i in range(nbatch)]}
chain_keys = ["nbatches", "ibatch"]
chain_configs = dict(
    nbatches,
    **ibatches,
) if nbatch > 1 else {}

# Set base directories to aggregate
run_groups = ['dt_rga']
channels   = ['ppim']
methods    = ["Asym", "HB"]
base_dirs  = [
    os.path.abspath(
        os.path.join(
                RGA_LAMBDA_ANALYSIS_HOME,
                f'jobs/saga/test_getKinBinned{method}__{rg}__{ch}__1D/'
            )
        ) for rg in run_groups for ch in channels for method in methods
]
YAML_DIR = os.path.abspath(
    os.path.join(
        RGA_LAMBDA_ANALYSIS_HOME,
        'yamls'
    )
)

# Set bin migration correction option
use_bin_mig = args.use_bin_mig

# Set list of channels for each base directory
chs = [ch for rg in run_groups for ch in channels for method in methods]

# Set channel label for each base directory
ch_sgasym_labels = {
    'ppim':'$D^{\\Lambda}_{LL\'}$',
}
ch_sgasym_labels = [ch_sgasym_labels[ch] for rg in run_groups for ch in channels for method in methods]

# Set x-axis labels for kinematic variables in all channels
xlabel_map = {
    'Q2':'$Q^{2}$ (GeV$^{2}$)', 'W':'$W$ (GeV)', 'y':'$y$', 'x':'$x$', 
    'z_ppim':'$z_{p\\pi^{-}}$', 'xF_ppim':'$x_{F p\\pi^{-}}$',
    'mass_ppim':'$M_{p\\pi^{-}}$ (GeV)',
}

# Set up list of run groups
rgs = [ rg for rg in run_groups for ch in channels for method in methods]

# Loop base directories
for rg, ch, base_dir, ch_sgasym_label in zip(rgs,chs,base_dirs,ch_sgasym_labels):

    # Setup input paths
    submit_path  = os.path.join(base_dir,"submit.sh")
    yaml_path    = os.path.join(base_dir,"args.yaml")
    out_path     = os.path.join(base_dir,"jobs.txt")
    #NOTE: Set the bin migration path below since this is binscheme dependent

    # Set aggregate keys
    aggregate_keys = []

    # # Load the binschemes from the path specified in the job yaml assuming there is only one given path and it is an absolute path
    # binschemes_paths_name = "binschemes_paths"
    # yaml_path = load_yaml(yaml_path)[binschemes_names][0]
    # binschemes = load_yaml(yaml_path)

    # Arguments for sagas.get_config_list()
    result_name = "a0" #NOTE: This also gets recycled as the asymmetry name

    # Arguments for sagas.get_out_dirs_list()
    sep='_'
    ext='.pdf'

    # Arguments for sagas.get_out_file_name()
    out_file_name_ext = '.csv'
    bin_mig_base_name="bin_mig_mat_"

    # Arguments for sagas.apply_bin_mig()
    # use_bin_mig = False
    id_gen_key='binid_gen'
    id_rec_key='binid_rec'
    mig_key='mig'
    results_keys = [result_name] #NOTE: You can apply bin migration to multiple dataframe entries in one go.

    # Arguments for sagas.get_graph_data()
    id_key = 'bin_id'

    # Arguments for sagas.get_graph_data()
    count_key  = 'count'
    asym_key   = result_name #NOTE: This is set from above
    err_ext    = '_err'

    # Arguments for saga.plot.plot_results()
    plot_results_kwargs_base = {
        'ylims':[-0.2,0.2],
        'sgasyms':[0.0],
        'sgasym_idx':0,
        'sgasym_labels':[ch_sgasym_label],
        'sg_colors':['blue'],
        'bgasyms':[],
        'bgasym_labels':[],
        'bg_colors':[],
        'show_injected_asymmetries':True,
        'hist_paths':[],
        'hist_colors':[
            'tab:red',
            'tab:blue',
        ],
        'hist_keys':[], #TODO: Set this and bin limits below...
        'hist_labels':[
            'RGA MC',
            'RGA Data',
        ],
        'watermark':'CLAS12 Preliminary',
        'hist_clone_axis':True,
        'hist_ylims':[0.0,0.06],
        'old_dat_path':None
    }

    # Additional useful parameters for plotting
    figsize = (16,10)
    #NOTE: Set outpath within the loop for unique naming
    use_default_plt_settings = True

    # If you want to rescale your results using results from other base directories set the following arguments
    rescale = False
    if rescale:
        plot_results_kwargs_base = dict(
            plot_results_kwargs_base,
            **{
                'old_dat_path':'',
                'new_sim_path':'',
                'old_sim_path':'',
                'count_key':'count',
                'yerr_key':'',
                'xs_ratio':1.0,
                'lumi_ratio':0.0,
            },
        )

    #---------- Set configurations ----------#
    # Setup configuration dictionary

    # Initialize configs
    configs = None
    aliases = None

    # # Load the binschemes from the path specified in the job yaml assuming there is only one given path and it is an absolute path
    # binschemes_paths_name = "binschemes_paths"
    # binscheme_yaml_path = load_yaml(yaml_path)[binschemes_paths_name][0]
    # binschemes = load_yaml(binscheme_yaml_path)
    binschemes = None

    # Aggregate with crystal ball signal
    if args.basic:

        # Set binschemes yaml path
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
        configs = dict(
            asymfitvars,
            **binschemes,
        )

        # Reset binschemes
        binschemes  = load_yaml(binscheme_yaml_path)

    # Aggregate with signal pdf types
    if args.massfit_types is not None:

        # Set binschemes yaml path
        binscheme_yaml_path = os.path.join(YAML_DIR,f'out_1d_bins_{ch}.yaml')

        # Create job submission structure
        binschemes  = load_yaml(binscheme_yaml_path)
        asymfitvars = {"asymfitvars":args.asymfitvars}
        binschemes  = {"binschemes":[{el:binschemes[el]} for el in binschemes]}

        # Create list of mass fit yaml file maps
        massfit_yamlfile_maps = [
            {
                f"scheme_{binscheme}_bin_{binid}": \
                os.path.join(
                    YAML_DIR,
                    f"massfit/{rg}/{massfit_type}/",
                    f"scheme_{binvar}_bin_{binid}.yaml",
                ) for binid in range(len(get_binscheme_cuts_and_ids(binscheme)[2]))
                for binscheme in binschemes 
            } for massfit_type in args.massfit_types
        ]

        # Print mass fit yaml file maps
        for idx, massfit_yamlfile_map in enumerate(massfit_yamlfile_maps):
            print("INFO: massfit_yamlfile_maps["+idx+"] = {")
            for key in massfit_yamlfile_map:
                print(f"INFO: \t{key}: {massfit_yamlfile_map[key]},")
            print("INFO: }")
            massfit_yamlfile_maps = {
                "massfit_yamlfile_map": massfit_yamlfile_maps
            }

        # Set aliases
        aliases     = {
            "binschemes":{
                str(el):list(el.keys())[0]+"_binscheme"
                for el in binschemes["binschemes"]
            },
            "massfit_yamlfile_map":{
                str(massfit_yamlfile_map):massfit_type \
                for massfit_type, massfit_yamlfile_map in \
                zip(args.massfit_types,massfit_yamlfile_maps)
            },
        }

        # Set replacements
        replacements = None

        # Set job file paths and configs
        configs = dict(
            asymfitvars,
            **binschemes,
            **massfit_yamlfile_maps,
        )

        # Reset bin schemes
        binschemes  = load_yaml(binscheme_yaml_path)

    # Aggregate with 2 cos phi regions
    if args.cos_phi:

        # Set binscheme yaml path
        binscheme_yaml_path = os.path.join(YAML_DIR,f'out_full_bin_{ch}.yaml')

        # Create job submission structure with fit variables and cos phi cuts
        asymfitvars = {"asymfitvars":args.asymfitvars}
        args_yaml_path = os.path.join(base_dir,"args.yaml")
        args_yaml = load_yaml(args_yaml_path)
        cuts = args_yaml["cuts"]
        cuts_pos_cos_phi = cuts + " && !(phi_h_ppim<TMath::Pi()/2 || phi_h_ppim>=3*TMath::Pi()/2)"
        cuts_neg_cos_phi = cuts + " && (phi_h_ppim<TMath::Pi()/2 || phi_h_ppim>=3*TMath::Pi()/2)"
        cutss = {"cuts":[cuts_pos_cos_phi,cuts_neg_cos_phi]}
        aliases     = {
            "cuts":{
                cuts_pos_cos_phi:"pos_cos_phi",
                cuts_neg_cos_phi:"neg_cos_phi",
            }
        }

        # Set job file paths and configs
        configs = dict(
            asymfitvars,
            **cutss,
        )

    # Aggregate with many bins and no sb subtraction binned in mass_ppim
    if args.no_sbs_mass_ppim:

        # Set binscheme yaml path
        binscheme_yaml_path = os.path.join(YAML_DIR,f"out_1d_bins_{ch}_no_bg_correction.yaml")

        # Create job submission structure with fit variables and mass_ppim binscheme and no sideband subtraction
        asymfitvars = {"asymfitvars":args.asymfitvars}
        use_sb_subtractions = {"use_sb_subtraction": [False]}
        # binscheme_yaml_path = os.path.join(YAML_DIR,f"out_1d_bins_{ch}_no_bg_correction.yaml")
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
        configs = dict(
            asymfitvars,
            **binschemes,
            **use_sb_subtractions,
        )

        # Reset binschemes
        binschemes = load_yaml(binscheme_yaml_path)["mass_ppim"]

    # Aggregate with many bins and no sb subtraction binned in xF_ppim with no xF>0 cut
    if args.no_sbs_xF_ppim:

        # Set binscheme yaml path
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
        configs = dict(
            asymfitvars,
            **binschemes,
            **use_sb_subtractions,
        )

        # Reset binschemes
        binschemes = load_yaml(binscheme_yaml_path)["xF_ppim"]

    # Aggregate with 1 bin and splot
    if args.splot:

        # Set binscheme yaml path
        binscheme_yaml_path = os.path.join(YAML_DIR,f'out_1d_bins_{ch}.yaml')

        # Create job submission structure with fit variables and cos phi cuts
        binschemes  = load_yaml(binscheme_yaml_path)
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
        configs = dict(
            asymfitvars,
            **binschemes,
            **use_splot,
        )

        # Reset binschemes
        binschemes  = load_yaml(binscheme_yaml_path)

    else:
        print("INFO: No usable configuration found, exiting.")
        sys.exit(0)

    # Check configs
    if configs is None:
        print("INFO: No usable configuration found, exiting.")
        sys.exit(0)

    # Get list of configurations
    config_list = sagas.get_config_list(configs,aggregate_keys=aggregate_keys)

    # Get aggregated list of directories
    out_dirs_list = sagas.get_out_dirs_list(
                                    configs,
                                    base_dir,
                                    aggregate_keys=aggregate_keys
                                )

    #---------- Loop bin schemes ----------#
    for binscheme_idx, binscheme_name in enumerate(binschemes.keys()):

        # Get the bin scheme
        binscheme = binschemes[binscheme_name]
        proj_var  = list(binscheme.keys())[0] #NOTE: Assume projection variable is the only variable in the bin scheme
        nbins = len(binscheme[proj_var])-1

        # Arguments for sagas.get_graph_data()
        xvar_keys = [proj_var]

        # Set some bin scheme dependent plotting parameters
        binlims = binscheme[proj_var]
        plot_results_kwargs_base['xlims'] = [binlims[0],binlims[-1]]
        plot_results_kwargs_base['xlabel'] = xlabel_map[binscheme_name]
        plot_results_kwargs_base['binlims'] = binlims
        plot_results_kwargs_base['hist_paths'] = [
            os.path.abspath(os.path.join(RGA_LAMBDA_ANALYSIS_HOME,f'jobs/saga/test_getBinKinematicsTH1Ds__mc_rga__{ch}/out_1d_bins_binscheme_kinematics.root')),
            os.path.abspath(os.path.join(RGA_LAMBDA_ANALYSIS_HOME,f'jobs/saga/test_getBinKinematicsTH1Ds__dt_rga__{ch}/out_1d_bins_binscheme_kinematics.root')),
        ]
        plot_results_kwargs_base['hist_keys'] = ['h1_bin0_'+binscheme_name for i in range(len(plot_results_kwargs_base['hist_paths']))]

        # Get the bin migration path
        bin_mig_base_dir = os.path.abspath(
            os.path.join(
                RGA_LAMBDA_ANALYSIS_HOME,
                f'jobs/saga/test_getBinMigration__{rg}__{ch}'
            )
        )
        bin_mig_path = sagas.get_out_file_name(
            base_dir=bin_mig_base_dir,
            base_name=bin_mig_base_name,
            binscheme_name=binscheme_name,
            ext=out_file_name_ext
        )

        # Load bin migration matrix and invert
        bin_mig_df, bin_mig_mat, inv_bin_mig_mat = None, None, None
        if use_bin_mig:
            bin_mig_df = load_csv(bin_mig_path)
            bin_mig_mat = sagas.get_bin_mig_mat(
                bin_mig_df,
                id_gen_key=id_gen_key,
                id_rec_key=id_rec_key,
                mig_key=mig_key,
            )
            save_bin_mig_mat_to_csv(
                bin_mig_mat,
                base_dir='./',
                basename=binscheme_name,
                delimiter=",",
                header=None,
                fmt=None,
                comments='',
            )
            inv_bin_mig_mat = np.linalg.inv(bin_mig_mat)

        #---------- Loop configurations ----------#
        # Loop each aggregate list
        for config_idx in range(len(config_list)):

            # Set the config you are interested in
            config = config_list[config_idx]
            out_dirs = out_dirs_list[config_idx]

            # Set the output path basename for this config
            config_out_path = sagas.get_config_out_path(
                    base_dir,
                    aggregate_keys,
                    binscheme_name+sep+rgh_mc_name+sep+result_name,
                    config,
                    sep=sep,
                    ext=ext,
                )
            config_out_path = os.path.join(base_dir,config_out_path)

            # Get the name of the CSV file for the binning scheme you are interested in
            out_file_names = [sagas.get_out_file_name(
                    base_dir=outdir,
                    base_name='out_',
                    binscheme_name=binscheme_name,
                    ext=out_file_name_ext
                ) for outdir in out_dirs]

            # Load pandas dataframes from the files
            dfs = [load_csv(out_file_name,config=config,chain_configs=chain_configs) for out_file_name in out_file_names]

            # Apply bin migration correction
            if use_bin_mig:
                for df in dfs:
                    sagas.apply_bin_mig(df,inv_bin_mig_mat,results_keys=results_keys) #NOTE: THIS MODIFIES THE DATAFRAMES IN PLACE

            # Get an aggregate graph
            proj_ids = [i for i in range(nbins)]#NOTE: Assume bin scheme indices are simple
            sgasym_idx = plot_results_kwargs_base['sgasym_idx'] #NOTE: Assume this is in the kwargs base dictionary
            sgasym = config['sgasyms'][sgasym_idx] if 'sgasyms' in config else 0.0
            plot_results_kwargs_base['sgasyms'] = config['sgasyms']#NOTE: TODO:
            aggregate_graph = sagas.get_aggregate_graph(
                [
                    sagas.get_graph_data(
                                df,
                                proj_ids,
                                id_key=id_key,
                                count_key=count_key,
                                xvar_keys=xvar_keys,
                                asym_key=asym_key,
                                err_ext=err_ext
                    ) for df in dfs
                ],
                xvar_keys=xvar_keys,
                sgasym=sgasym
            )

            # Use default plotting settings
            if use_default_plt_settings: set_default_plt_settings()

            # Create figure and axes
            f, ax = plt.subplots(figsize=figsize)

            # Set additional arguments for saga.plot.plot_results()
            plot_results_kwargs_base['sgasyms'] = config['sgasyms']
            plot_results_kwargs_base['outpath'] = config_out_path

            # Plot the graph
            plot_results(ax,**aggregate_graph,**plot_results_kwargs_base)

            # Save the graph
            f.savefig(config_out_path)

            # Close the graph
            plt.close()
